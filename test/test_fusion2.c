/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#ifndef P4_TO_P8
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_bits.h>
#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 4;
#endif

typedef struct
{
  int                 dummy;
}
test_fusion_t;

enum
{
  TIMINGS_INIT,
  TIMINGS_NUM_STATS
};

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;
  int                 eps = 0.05 * (P4EST_ROOT_LEN / 2);
  int                 center = (quadrant->x - (P4EST_ROOT_LEN / 2) ^ 2 +
                                quadrant->y - (P4EST_ROOT_LEN / 2) ^ 2 +
#ifdef P4_TO_P8
                                quadrant->z - (P4EST_ROOT_LEN / 2) ^ 2
#endif
    );

  /* Don't refine deeper than a given maximum level. */
  if (quadrant->level > refine_level) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  /* Trying to emulate refining edge around the circle *
   * with margin being 0.05 of the ROOT_LEN.               */
  if (center <= (P4EST_ROOT_LEN / 2 + eps) ^ 2 &&
      center >= (P4EST_ROOT_LEN / 2 - eps) ^ 2) {
    return 1;
  }

  /* Filter out large enough quadrants to refine.
   * && or || ??                                        */
  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - refine_level) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - refine_level)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - refine_level)
#endif
      )) {
    return 1;
  }

  /* from ghost. what is which_tree ?? */
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }

  return 1;
}

typedef struct
{
  int                 counter;
  int                *refine_flags;
}
refine_loop_t;

/* create an array of ints with values meaning:
 *   0 : keep
 *   1 : refine
 *   2 : coarsen
 */
enum
{
  FUSION_KEEP = 0,
  FUSION_REFINE,
  FUSION_COARSEN
};

static int
refine_in_loop (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quad)
{
  int                 flag;
  refine_loop_t      *loop_ctx = (refine_loop_t *) (p4est->user_pointer);

  flag = (loop_ctx->refine_flags[loop_ctx->counter++]);
  if (flag == FUSION_REFINE) {
    return 1;
  }
  return 0;
}

static void
mark_leaves (p4est_t * p4est, int *refine_flags, test_fusion_t * ctx)
{
  p4est_locidx_t      i;
  p4est_locidx_t      num_local = p4est->local_num_quadrants;

  /* TODO: replace with a meaningful (or more than one meaningful, controlled
   * by ctx) refinment pattern */
  for (i = 0; i < num_local; i++) {
    /* refine every other */
    refine_flags[i] = (i & 1) ? FUSION_REFINE : FUSION_KEEP;
  }
}

enum
{
  FUSION_FULL_LOOP,
  FUSION_TIME_COARSEN,
  FUSION_TIME_REFINE,
  FUSION_TIME_BALANCE,
  FUSION_TIME_PARTITION,
  FUSION_TIME_GHOST,
  FUSION_NUM_STATS
};

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  int                 i;
  p4est_lnodes_t     *lnodes;   /* runtime option: time lnodes construction */
  int                 first_argc;
  int                 type;
  int                 num_tests = 3;
  sc_statinfo_t       stats[FUSION_NUM_STATS];
  sc_options_t       *opt;
  int                 log_priority = SC_LP_ESSENTIAL;
  test_fusion_t       ctx;

  /* initialize MPI and libsc, p4est packages */
  mpiret = sc_MPI_Init (&argc, &argv);
  mpicomm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'k', "num-tests", &num_tests, num_tests,
                      "The number of instances of timiing the fusion loop");
  sc_options_add_int (opt, 'q', "quiet", &log_priority, log_priority,
                      "Degree of quietude of output");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  sc_set_log_defaults (NULL, NULL, log_priority);
  p4est_init (NULL, log_priority);

  /* set values in ctx here */
  /* random 1 for now */
  ctx.dummy = 1;

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  p4est_refine (p4est, 1, refine_fn, NULL);

  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  p4est_partition (p4est, 0, NULL);

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  sc_stats_init (&stats[FUSION_FULL_LOOP], "Full loop");
  sc_stats_init (&stats[FUSION_TIME_REFINE], "Local refinement");
  sc_stats_init (&stats[FUSION_TIME_BALANCE], "Balance");
  sc_stats_init (&stats[FUSION_TIME_PARTITION], "Partition");
  sc_stats_init (&stats[FUSION_TIME_GHOST], "Ghost");

  for (i = 0; i <= num_tests; i++) {
    p4est_t            *forest_copy;
    int                *refine_flags, *rflags_copy;
    p4est_ghost_t      *gl_copy;
    refine_loop_t       loop_ctx;
    sc_flopinfo_t       fi_full, snapshot_full;
    sc_flopinfo_t       fi_refine, snapshot_refine;
    sc_flopinfo_t       fi_balance, snapshot_balance;
    sc_flopinfo_t       fi_partition, snapshot_partition;
    sc_flopinfo_t       fi_ghost, snapshot_ghost;
    sc_flopinfo_t       fi_coarsen, snapshot_coarsen;

    if (!i) {
      P4EST_GLOBAL_PRODUCTION ("Timing loop 0 (discarded)\n");
    }
    else {
      P4EST_GLOBAL_PRODUCTIONF ("Timing loop %d\n", i);
    }
    sc_log_indent_push_count (p4est_package_id, 2);

    forest_copy = p4est_copy (p4est, 0 /* do not copy data */ );

    /* predefine which leaves we want to refine and coarsen */

    refine_flags = P4EST_ALLOC (int, p4est->local_num_quadrants);

    mark_leaves (forest_copy, refine_flags, &ctx);

    /* Once the leaves have been marked, we have sufficient information to
     * complete a refinement cycle: start the timing at this point */
    sc_flops_snap (&fi_full, &snapshot_full);

    /* start the timing of one instance of the timing cycle */
    /* see sc_flops_snap() / sc_flops_shot() in timings2.c */

    /* non-recursive refinement loop: the callback simply checks the flags
     * that we have defined for which leaves we want to refine */

    /* TODO: conduct coarsening before refinement.  This requires making a
     * copy of refine_flags: creating one that is valid after coarsening */

    sc_flops_snap (&fi_coarsen, &snapshot_coarsen);

    /* coarsen */
    rflags_copy = P4EST_STRDUP (refine_flags);

    sc_flops_shot (&fi_coarsen, &snapshot_coarsen);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_COARSEN],
                           snapshot_coarsen.iwtime);
    }

    sc_flops_snap (&fi_refine, &snapshot_refine);
    loop_ctx.counter = 0;
    loop_ctx.refine_flags = refine_flags;

    forest_copy->user_pointer = (void *) &loop_ctx;
    p4est_refine (forest_copy, 0 /* non-recursive */ , refine_in_loop, NULL);
    sc_flops_shot (&fi_refine, &snapshot_refine);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_REFINE],
                           snapshot_refine.iwtime);
    }

    sc_flops_shot (&fi_balance, &snapshot_balance);
    p4est_balance (forest_copy, P4EST_CONNECT_FULL, NULL);
    sc_flops_shot (&fi_balance, &snapshot_balance);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_BALANCE],
                           snapshot_balance.iwtime);
    }

    sc_flops_shot (&fi_partition, &snapshot_partition);
    p4est_partition (forest_copy, 0, NULL);
    sc_flops_shot (&fi_partition, &snapshot_partition);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_PARTITION],
                           snapshot_partition.iwtime);
    }

    sc_flops_shot (&fi_ghost, &snapshot_ghost);
    gl_copy = p4est_ghost_new (forest_copy, P4EST_CONNECT_FULL);
    sc_flops_shot (&fi_ghost, &snapshot_ghost);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_GHOST], snapshot_ghost.iwtime);
    }

    /* end  the timing of one instance of the timing cycle */
    sc_flops_shot (&fi_full, &snapshot_full);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_FULL_LOOP], snapshot_full.iwtime);
    }

    /* clean up */
    P4EST_FREE (refine_flags);
    p4est_ghost_destroy (gl_copy);
    p4est_destroy (forest_copy);

    sc_log_indent_pop_count (p4est_package_id, 2);
  }

  /* accumulate and print statistics */
  sc_stats_compute (mpicomm, FUSION_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  FUSION_NUM_STATS, stats, 1, 1);

  /* clean up */
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_options_destroy (opt);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
