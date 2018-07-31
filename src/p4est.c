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

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_io.h>
#include <p8est_search.h>
#include <p8est_balance.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_io.h>
#include <p4est_search.h>
#include <p4est_balance.h>
#endif /* !P4_TO_P8 */
#include <sc_io.h>
#include <sc_notify.h>
#include <sc_ranges.h>
#include <sc_search.h>
#ifdef P4EST_HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef P4EST_ENABLE_MPIIO
#define P4EST_MPIIO_WRITE
#endif

#ifdef P4EST_HAVE_UNISTD_H
#include <unistd.h>
#endif

typedef struct
{
  int8_t              have_first_count, have_first_load;
  int8_t              have_second_count, have_second_load;
  int                 recv_first_count, recv_second_count;
  int                 send_first_count, send_second_count;
  sc_array_t          send_first, send_second, recv_first, recv_second;
}
p4est_balance_peer_t;

#define p4est_num_ranges (25)

#ifndef P4_TO_P8

static int          p4est_uninitialized_key;
void               *P4EST_DATA_UNINITIALIZED = &p4est_uninitialized_key;

#endif /* P4_TO_P8 */

static const size_t number_toread_quadrants = 32;
static const int8_t fully_owned_flag = 0x01;
static const int8_t any_face_flag = 0x02;

void
p4est_qcoord_to_vertex (p4est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
#ifdef P4_TO_P8
                        p4est_qcoord_t z,
#endif
                        double vxyz[3])
{
  const double       *vertices = connectivity->vertices;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_topidx_t num_vertices = connectivity->num_vertices;
#endif
  const p4est_topidx_t *vindices;
  int                 xi, yi;
  double              wx[2], wy[2];
#ifdef P4_TO_P8
  int                 zi;
  double              wz[2];
#endif
  double              xfactor, yfactor;
  p4est_topidx_t      vindex;

  P4EST_ASSERT (num_vertices > 0);
  P4EST_ASSERT (vertices != NULL);
  P4EST_ASSERT (treeid >= 0 && treeid < connectivity->num_trees);

  P4EST_ASSERT (connectivity->tree_to_vertex != NULL);
  vindices = connectivity->tree_to_vertex + P4EST_CHILDREN * treeid;

  vxyz[0] = vxyz[1] = vxyz[2] = 0.;

  P4EST_ASSERT (x >= 0 && x <= P4EST_ROOT_LEN);
  wx[1] = (double) x / (double) P4EST_ROOT_LEN;
  wx[0] = 1. - wx[1];

  P4EST_ASSERT (y >= 0 && y <= P4EST_ROOT_LEN);
  wy[1] = (double) y / (double) P4EST_ROOT_LEN;
  wy[0] = 1. - wy[1];

#ifdef P4_TO_P8
  P4EST_ASSERT (z >= 0 && z <= P4EST_ROOT_LEN);
  wz[1] = (double) z / (double) P4EST_ROOT_LEN;
  wz[0] = 1. - wz[1];

  for (zi = 0; zi < 2; ++zi) {
#endif
    for (yi = 0; yi < 2; ++yi) {
#ifdef P4_TO_P8
      yfactor = wz[zi] * wy[yi];
#else
      yfactor = wy[yi];
#endif
      for (xi = 0; xi < 2; ++xi) {
        xfactor = yfactor * wx[xi];

        vindex = *vindices++;
        P4EST_ASSERT (vindex >= 0 && vindex < num_vertices);

        vxyz[0] += xfactor * vertices[3 * vindex + 0];
        vxyz[1] += xfactor * vertices[3 * vindex + 1];
        vxyz[2] += xfactor * vertices[3 * vindex + 2];
      }
    }
#ifdef P4_TO_P8
  }
#endif
}

size_t
p4est_memory_used (p4est_t * p4est)
{
  int                 mpisize;
  size_t              size;
  p4est_topidx_t      nt;
  p4est_tree_t       *tree;

  /* do not assert p4est_is_valid since it is collective */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);
  P4EST_ASSERT (p4est->trees != NULL);

  mpisize = p4est->mpisize;
  size = sizeof (p4est_t) +
    (mpisize + 1) * (sizeof (p4est_gloidx_t) + sizeof (p4est_quadrant_t));

  size += sc_array_memory_used (p4est->trees, 1);
  for (nt = 0; nt < p4est->connectivity->num_trees; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    size += sc_array_memory_used (&tree->quadrants, 0);
  }

  if (p4est->data_size > 0) {
    P4EST_ASSERT (p4est->user_data_pool != NULL);
    size += sc_mempool_memory_used (p4est->user_data_pool);
  }
  P4EST_ASSERT (p4est->quadrant_pool != NULL);
  size += sc_mempool_memory_used (p4est->quadrant_pool);

  return size;
}

long
p4est_revision (p4est_t * p4est)
{
  /* do not assert p4est_is_valid since it is collective */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->revision >= 0);

  return p4est->revision;
}

p4est_t            *
p4est_new (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
           size_t data_size, p4est_init_t init_fn, void *user_pointer)
{
  return p4est_new_ext (mpicomm, connectivity, 0, 0, 1,
                        data_size, init_fn, user_pointer);
}

p4est_t            *
p4est_new_ext (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
               p4est_locidx_t min_quadrants, int min_level, int fill_uniform,
               size_t data_size, p4est_init_t init_fn, void *user_pointer)
{
  int                 num_procs, rank;
  int                 i, must_remove_last_quadrant;
  int                 level;
  uint64_t            first_morton, last_morton, miu, count;
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t      tree_num_quadrants, global_num_quadrants;
  p4est_gloidx_t      first_tree, first_quadrant, first_tree_quadrant;
  p4est_gloidx_t      last_tree, last_quadrant, last_tree_quadrant;
  p4est_gloidx_t      quadrant_index;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_quadrant_t    a, b, c;
  p4est_quadrant_t   *global_first_position;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into " P4EST_STRING
     "_new with min quadrants %lld level %d uniform %d\n",
     (long long) min_quadrants, SC_MAX (min_level, 0), fill_uniform);
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));
  P4EST_ASSERT (min_level <= P4EST_QMAXLEVEL);

  /* create p4est object and assign some data members */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  p4est->data_size = data_size;
  p4est->user_pointer = user_pointer;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* set parallel environment */
  p4est_comm_parallel_env_assign (p4est, mpicomm);
  num_procs = p4est->mpisize;
  rank = p4est->mpirank;

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    if (tree_num_quadrants >=
        (num_procs * (p4est_gloidx_t) min_quadrants + (num_trees - 1))
        / num_trees) {
      break;
    }
    tree_num_quadrants *= P4EST_CHILDREN;
    P4EST_ASSERT (tree_num_quadrants > 0);
  }
  for (; level < min_level; ++level) {
    tree_num_quadrants *= P4EST_CHILDREN;
    P4EST_ASSERT (tree_num_quadrants > 0);
  }
  P4EST_ASSERT (level <= P4EST_QMAXLEVEL
                && tree_num_quadrants <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);

  /* compute global number of quadrants */
  global_num_quadrants = tree_num_quadrants * num_trees;
  P4EST_GLOBAL_PRODUCTIONF ("New " P4EST_STRING
                            " with %lld trees on %d processors\n",
                            (long long) num_trees, num_procs);
  P4EST_GLOBAL_INFOF ("Initial level %d potential global quadrants"
                      " %lld per tree %lld\n",
                      level, (long long) global_num_quadrants,
                      (long long) tree_num_quadrants);

  /* compute index of first tree for this processor */
  first_quadrant =
    p4est_partition_cut_gloidx (global_num_quadrants, rank, num_procs);
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant = p4est_partition_cut_gloidx (global_num_quadrants, rank + 1,
                                              num_procs) - 1;
  P4EST_VERBOSEF
    ("first tree %lld first quadrant %lld global quadrant %lld\n",
     (long long) first_tree, (long long) first_tree_quadrant,
     (long long) first_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);

  /* compute index of last tree for this processor */
  if (first_quadrant <= last_quadrant) {
    last_tree = last_quadrant / tree_num_quadrants;
    last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;
    P4EST_VERBOSEF
      ("last tree %lld last quadrant %lld global quadrant %lld\n",
       (long long) last_tree, (long long) last_tree_quadrant,
       (long long) last_quadrant);

    /* check ranges of various integers to be 32bit compatible */
    P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
    P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
    P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);
    if (first_tree == last_tree) {
      P4EST_ASSERT (first_tree_quadrant <= last_tree_quadrant);
    }
  }
  else {
    P4EST_VERBOSE ("Empty processor");
    P4EST_ASSERT (0 <= first_tree && 0 <= first_tree_quadrant);
    first_tree = -1;
    last_tree = -2;
    last_tree_quadrant = -1;
  }

  /* allocate trees and quadrants */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
    tree->quadrants_offset = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = -1;
    }
    tree->maxlevel = 0;
  }
  p4est->local_num_quadrants = 0;
  p4est->global_num_quadrants = 0;

  /* for every locally non-empty tree fill first and last quadrant */
  P4EST_QUADRANT_INIT (&a);
  P4EST_QUADRANT_INIT (&b);
  P4EST_QUADRANT_INIT (&c);
  for (jt = first_tree; jt <= last_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;

    quad = NULL;
    if (!fill_uniform) {        /* fill with coarsest possible quadrants */
      must_remove_last_quadrant = 0;

      /* set morton id of first quadrant and initialize user data */
      if (jt == first_tree) {
        p4est_quadrant_set_morton (&a, level, first_tree_quadrant);
      }
      else {
        p4est_quadrant_set_morton (&a, level, 0);
      }
#ifdef P4_TO_P8
      P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx 0x%llx\n",
                     (long long) jt, (long long) a.x,
                     (long long) a.y, (long long) a.z);
#else
      P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx\n",
                     (long long) jt, (long long) a.x, (long long) a.y);
#endif
      p4est_quadrant_first_descendant (&a, &tree->first_desc,
                                       P4EST_QMAXLEVEL);

      /* set morton id of last quadrant */
      if (tree_num_quadrants == 1 ||
          (jt == first_tree
           && first_tree_quadrant == tree_num_quadrants - 1)) {
        /* There is only a in the tree */
        quad = p4est_quadrant_array_push (tquadrants);
        *quad = a;
        p4est_quadrant_init_data (p4est, jt, quad, init_fn);
        tree->maxlevel = a.level;
        tree->quadrants_per_level[a.level] = 1;
      }
      else {
        if (jt == last_tree) {
          if (last_tree_quadrant == tree_num_quadrants - 1) {
            quadrant_index = last_tree_quadrant;
          }
          else {
            quadrant_index = last_tree_quadrant + 1;
            must_remove_last_quadrant = 1;
          }
          p4est_quadrant_set_morton (&b, level, quadrant_index);
        }
        else {
          p4est_quadrant_set_morton (&b, level, tree_num_quadrants - 1);
        }
#ifdef P4_TO_P8
        P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx 0x%llx\n",
                       (long long) jt, (long long) b.x,
                       (long long) b.y, (long long) b.z);
#else
        P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx\n",
                       (long long) jt, (long long) b.x, (long long) b.y);
#endif
        /* fill up tree between a and b with coarse quadrants */
        p4est_complete_region (p4est, &a, 1, &b, !must_remove_last_quadrant,
                               tree, jt, init_fn);
        quad = p4est_quadrant_array_index (tquadrants,
                                           tquadrants->elem_count - 1);
      }
    }
    else {                      /* fill tree with quadrants of given level */
      /* determine range of quadrants in this tree */
      first_morton = (uint64_t)
        (jt == first_tree ? first_tree_quadrant : 0);
      last_morton = (uint64_t)
        (jt == last_tree ? last_tree_quadrant : tree_num_quadrants - 1);
      count = last_morton - first_morton + 1;
      P4EST_ASSERT (count > 0);

      /* populate quadrant array in Morton order */
      sc_array_resize (tquadrants, (size_t) count);
      for (miu = 0; miu < count; ++miu) {
        quad = p4est_quadrant_array_index (tquadrants, (size_t) miu);
        p4est_quadrant_set_morton (quad, level, first_morton + miu);
        p4est_quadrant_init_data (p4est, jt, quad, init_fn);
      }

      /* remember first tree position */
      p4est_quadrant_first_descendant (p4est_quadrant_array_index
                                       (tquadrants, 0), &tree->first_desc,
                                       P4EST_QMAXLEVEL);

      /* set tree counters */
      tree->maxlevel = (int8_t) level;
      tree->quadrants_per_level[level] = (p4est_locidx_t) count;
    }

#if 0
    P4EST_VERBOSEF ("tree %lld quadrants %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);
#endif

    tree->quadrants_offset = p4est->local_num_quadrants;
    p4est->local_num_quadrants += tquadrants->elem_count;
    p4est_quadrant_last_descendant (quad, &tree->last_desc, P4EST_QMAXLEVEL);
  }
  if (last_tree >= 0) {
    for (; jt < num_trees; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  /* compute some member variables */
  p4est->first_local_tree = first_tree;
  p4est->last_local_tree = last_tree;
  p4est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  if (!fill_uniform && level > 0) {
    /* this performs an allgather to count all quadrants */
    p4est_comm_count_quadrants (p4est);
  }
  else {
    /* for a uniform forest we know all global information a priori */
    for (i = 0; i <= num_procs; ++i) {
      p4est->global_first_quadrant[i] =
        p4est_partition_cut_gloidx (global_num_quadrants, i, num_procs);
    }
    p4est->global_num_quadrants = global_num_quadrants;
  }

  /* fill in global partition information */
  global_first_position = P4EST_ALLOC_ZERO (p4est_quadrant_t, num_procs + 1);
  for (i = 0; i <= num_procs; ++i) {
    first_quadrant =
      p4est_partition_cut_gloidx (global_num_quadrants, i, num_procs);
    first_tree = first_quadrant / tree_num_quadrants;
    first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
    p4est_quadrant_set_morton (&c, level, first_tree_quadrant);
    global_first_position[i].x = c.x;
    global_first_position[i].y = c.y;
#ifdef P4_TO_P8
    global_first_position[i].z = c.z;
#endif
    global_first_position[i].level = P4EST_QMAXLEVEL;
    global_first_position[i].p.which_tree = first_tree;
  }
  p4est->global_first_position = global_first_position;

  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est->revision == 0);
  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_new with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              qz;
#endif
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;

  for (jt = 0; jt < p4est->connectivity->num_trees; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);

#ifdef P4EST_ENABLE_DEBUG
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz) {
      p4est_quadrant_t   *quad =
        p4est_quadrant_array_index (&tree->quadrants, qz);
      p4est_quadrant_free_data (p4est, quad);
    }
#endif

    sc_array_reset (&tree->quadrants);
  }
  sc_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    sc_mempool_destroy (p4est->user_data_pool);
  }
  sc_mempool_destroy (p4est->quadrant_pool);

  p4est_comm_parallel_env_release (p4est);
  P4EST_FREE (p4est->global_first_quadrant);
  P4EST_FREE (p4est->global_first_position);
  P4EST_FREE (p4est);
}

p4est_t            *
p4est_copy (p4est_t * input, int copy_data)
{
  return p4est_copy_ext (input, copy_data, 0 /* don't duplicate MPI comm */ );
}

p4est_t            *
p4est_copy_ext (p4est_t * input, int copy_data, int duplicate_mpicomm)
{
  const p4est_topidx_t num_trees = input->connectivity->num_trees;
  const p4est_topidx_t first_tree = input->first_local_tree;
  const p4est_topidx_t last_tree = input->last_local_tree;
  size_t              icount;
  size_t              zz;
  p4est_topidx_t      jt;
  p4est_t            *p4est;
  p4est_tree_t       *itree, *ptree;
  p4est_quadrant_t   *iq, *pq;
  sc_array_t         *iquadrants, *pquadrants;

  /* create a shallow copy and zero out dependent fields */
  p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, input, sizeof (p4est_t));
  p4est->global_first_quadrant = NULL;
  p4est->global_first_position = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;

  /* set parallel environment */
  p4est_comm_parallel_env_assign (p4est, input->mpicomm);
  if (duplicate_mpicomm) {
    p4est_comm_parallel_env_duplicate (p4est);
  }

  /* allocate a user data pool if necessary and a quadrant pool */
  if (copy_data && p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->data_size = 0;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* copy quadrants for each tree */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    ptree = p4est_tree_array_index (p4est->trees, jt);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
  }
  for (jt = first_tree; jt <= last_tree; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    iquadrants = &itree->quadrants;
    icount = iquadrants->elem_count;
    ptree = p4est_tree_array_index (p4est->trees, jt);
    pquadrants = &ptree->quadrants;
    sc_array_resize (pquadrants, icount);
    memcpy (pquadrants->array, iquadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      P4EST_ASSERT (copy_data);
      for (zz = 0; zz < icount; ++zz) {
        iq = p4est_quadrant_array_index (iquadrants, zz);
        pq = p4est_quadrant_array_index (pquadrants, zz);
        pq->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        memcpy (pq->p.user_data, iq->p.user_data, p4est->data_size);
      }
    }
  }

  /* allocate and copy global quadrant count */
  p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize + 1);
  memcpy (p4est->global_first_quadrant, input->global_first_quadrant,
          (p4est->mpisize + 1) * sizeof (p4est_gloidx_t));

  /* allocate and copy global partition information */
  p4est->global_first_position = P4EST_ALLOC (p4est_quadrant_t,
                                              p4est->mpisize + 1);
  memcpy (p4est->global_first_position, input->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_quadrant_t));

  /* the copy starts with a revision count of zero */
  p4est->revision = 0;

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}

void
p4est_reset_data (p4est_t * p4est, size_t data_size,
                  p4est_init_t init_fn, void *user_pointer)
{
  int                 doresize;
  size_t              zz;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;

  doresize = (p4est->data_size != data_size);

  p4est->data_size = data_size;
  p4est->user_pointer = user_pointer;

  if (doresize) {
    if (p4est->user_data_pool != NULL) {
      sc_mempool_destroy (p4est->user_data_pool);
    }
    if (p4est->data_size > 0) {
      p4est->user_data_pool = sc_mempool_new (p4est->data_size);
    }
    else {
      p4est->user_data_pool = NULL;
    }
  }

  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      q = p4est_quadrant_array_index (tquadrants, zz);
      if (doresize) {
        if (p4est->data_size > 0) {
          q->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        }
        else {
          q->p.user_data = NULL;
        }
      }
      if (init_fn != NULL) {
        init_fn (p4est, jt, q);
      }
    }
  }
}

void
p4est_refine (p4est_t * p4est, int refine_recursive,
              p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  p4est_refine_ext (p4est, refine_recursive, -1, refine_fn, init_fn, NULL);
}

void
p4est_refine_ext (p4est_t * p4est, int refine_recursive, int allowed_level,
                  p4est_refine_t refine_fn, p4est_init_t init_fn,
                  p4est_replace_t replace_fn)
{
  p4est_gloidx_t      old_gnq;

  if (allowed_level < 0) {
    allowed_level = P4EST_QMAXLEVEL;
  }
  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_refine with %lld total quadrants,"
                            " allowed level %d\n",
                            (long long) p4est->global_num_quadrants,
                            allowed_level);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not decrease */
  old_gnq = p4est->global_num_quadrants;
  p4est_refine_ext_dirty (p4est, refine_recursive, allowed_level, refine_fn,
                          init_fn, replace_fn);
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  P4EST_ASSERT (p4est->global_num_quadrants >= old_gnq);
  if (old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_refine_ext_dirty (p4est_t * p4est, int refine_recursive,
                        int allowed_level, p4est_refine_t refine_fn,
                        p4est_init_t init_fn, p4est_replace_t replace_fn)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              quadrant_pool_size, data_pool_size;
#endif
  int                 firsttime;
  int                 i, maxlevel;
  p4est_topidx_t      nt;
  size_t              incount, current, restpos, movecount;
  sc_list_t          *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *family[8];
  sc_flopinfo_t       snap;
  p4est_quadrant_t    parent, *pp = &parent;

  P4EST_FUNC_SNAP (p4est, &snap);

  if (allowed_level < 0) {
    allowed_level = P4EST_QMAXLEVEL;
  }
  P4EST_ASSERT (0 <= allowed_level && allowed_level <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (refine_fn != NULL);

  /*
     q points to a quadrant that is an array member
     qalloc is a quadrant that has been allocated through quadrant_pool
     qpop is a quadrant that has been allocated through quadrant_pool
     never mix these two types of quadrant pointers

     The quadrant->pad8 field of list quadrants is interpreted as boolean
     and set to true for quadrants that have already been refined.
   */
  list = sc_list_new (NULL);
  p4est->local_num_quadrants = 0;

  /* loop over all local trees */
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    tree->quadrants_offset = p4est->local_num_quadrants;
    tquadrants = &tree->quadrants;
#ifdef P4EST_ENABLE_DEBUG
    quadrant_pool_size = p4est->quadrant_pool->elem_count;
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
#endif

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into refine tree %lld with %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);

    /* reset the quadrant counters */
    maxlevel = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    incount = tquadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = p4est_quadrant_array_index (tquadrants, current);
      if (refine_fn (p4est, nt, q) && (int) q->level < allowed_level) {
        break;
      }
      maxlevel = SC_MAX (maxlevel, (int) q->level);
      ++tree->quadrants_per_level[q->level];
    }
    if (current == incount) {
      /* no refinement occurs in this tree */
      p4est->local_num_quadrants += incount;
      continue;
    }
    P4EST_ASSERT (q != NULL);

    /* now we have a quadrant to refine, prepend it to the list */
    qalloc = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
    *qalloc = *q;               /* never prepend array members directly */
    qalloc->pad8 = 0;           /* this quadrant has not been refined yet */
    (void) sc_list_prepend (list, qalloc);      /* only new quadrants */

    P4EST_QUADRANT_INIT (&parent);

    /*
       current points to the next array member to write
       restpos points to the next array member to read
     */
    restpos = current + 1;

    /* run through the list and refine recursively */
    firsttime = 1;
    while (list->elem_count > 0) {
      qpop = p4est_quadrant_list_pop (list);
      if (firsttime ||
          ((refine_recursive || !qpop->pad8) &&
           refine_fn (p4est, nt, qpop) &&
           (int) qpop->level < allowed_level)) {
        firsttime = 0;
        sc_array_resize (tquadrants,
                         tquadrants->elem_count + P4EST_CHILDREN - 1);

        if (replace_fn != NULL) {
          /* do not free qpop's data yet: we will do this when the parent
           * is replaced */
          parent = *qpop;
        }
        else {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c2 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c3 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);

#ifdef P4_TO_P8
        c4 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c5 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c6 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c7 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);

        p8est_quadrant_children (qpop, c0, c1, c2, c3, c4, c5, c6, c7);
#else
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
#endif
        p4est_quadrant_init_data (p4est, nt, c0, init_fn);
        p4est_quadrant_init_data (p4est, nt, c1, init_fn);
        p4est_quadrant_init_data (p4est, nt, c2, init_fn);
        p4est_quadrant_init_data (p4est, nt, c3, init_fn);
        c0->pad8 = c1->pad8 = c2->pad8 = c3->pad8 = 1;

#ifdef P4_TO_P8
        p4est_quadrant_init_data (p4est, nt, c4, init_fn);
        p4est_quadrant_init_data (p4est, nt, c5, init_fn);
        p4est_quadrant_init_data (p4est, nt, c6, init_fn);
        p4est_quadrant_init_data (p4est, nt, c7, init_fn);
        c4->pad8 = c5->pad8 = c6->pad8 = c7->pad8 = 1;

        (void) sc_list_prepend (list, c7);
        (void) sc_list_prepend (list, c6);
        (void) sc_list_prepend (list, c5);
        (void) sc_list_prepend (list, c4);
#endif
        (void) sc_list_prepend (list, c3);
        (void) sc_list_prepend (list, c2);
        (void) sc_list_prepend (list, c1);
        (void) sc_list_prepend (list, c0);

        if (replace_fn != NULL) {
          /* in family mode we always call the replace callback right
           * away */
          family[0] = c0;
          family[1] = c1;
          family[2] = c2;
          family[3] = c3;
#ifdef P4_TO_P8
          family[4] = c4;
          family[5] = c5;
          family[6] = c6;
          family[7] = c7;
#endif
          replace_fn (p4est, nt, 1, &pp, P4EST_CHILDREN, family);
          p4est_quadrant_free_data (p4est, &parent);
        }
      }
      else {
        /* need to make room in the array to store this new quadrant */
        if (restpos < incount && current == restpos) {
          movecount = SC_MIN (incount - restpos, number_toread_quadrants);
          while (movecount > 0) {
            q = p4est_quadrant_array_index (tquadrants, restpos);
            qalloc = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
            *qalloc = *q;       /* never append array members directly */
            qalloc->pad8 = 0;   /* has not been refined yet */
            (void) sc_list_append (list, qalloc);       /* only new quadrants */
            --movecount;
            ++restpos;
          }
        }

        /* store new quadrant and update counters */
        q = p4est_quadrant_array_index (tquadrants, current);
        *q = *qpop;
        maxlevel = SC_MAX (maxlevel, (int) qpop->level);
        ++tree->quadrants_per_level[qpop->level];
        ++current;
        sc_mempool_free (p4est->quadrant_pool, qpop);
      }
    }
    tree->maxlevel = (int8_t) maxlevel;
    p4est->local_num_quadrants += tquadrants->elem_count;

    P4EST_ASSERT (restpos == incount);
    P4EST_ASSERT (current == tquadrants->elem_count);
    P4EST_ASSERT (list->first == NULL && list->last == NULL);
    P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size + tquadrants->elem_count ==
                    p4est->user_data_pool->elem_count + incount);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done refine tree %lld now %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);
  }
  if (p4est->last_local_tree >= 0) {
    for (; nt < p4est->connectivity->num_trees; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  sc_list_destroy (list);

  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_coarsen (p4est_t * p4est, int coarsen_recursive,
               p4est_coarsen_t coarsen_fn, p4est_init_t init_fn)
{
  p4est_coarsen_ext (p4est, coarsen_recursive, 0, coarsen_fn, init_fn, NULL);
}

void
p4est_coarsen_ext (p4est_t * p4est, int coarsen_recursive,
                   int callback_orphans, p4est_coarsen_t coarsen_fn,
                   p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_gloidx_t      old_gnq;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not increase */
  old_gnq = p4est->global_num_quadrants;
  p4est_coarsen_ext_dirty (p4est, coarsen_recursive, callback_orphans,
                           coarsen_fn, init_fn, replace_fn);
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  P4EST_ASSERT (p4est->global_num_quadrants <= old_gnq);
  if (old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_coarsen_ext_dirty (p4est_t * p4est,
                         int coarsen_recursive, int callback_orphans,
                         p4est_coarsen_t coarsen_fn, p4est_init_t init_fn,
                         p4est_replace_t replace_fn)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  int                 i, maxlevel;
  int                 isfamily;
  size_t              zz;
  size_t              incount, removed;
  size_t              window, start, length, cidz;
  p4est_locidx_t      num_quadrants, prev_offset;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *c[P4EST_CHILDREN];
  p4est_quadrant_t   *cfirst, *clast;
  sc_array_t         *tquadrants;
  sc_flopinfo_t       snap;
  p4est_quadrant_t    qtemp;

  P4EST_FUNC_SNAP (p4est, &snap);

  P4EST_ASSERT (coarsen_fn != NULL);

  P4EST_QUADRANT_INIT (&qtemp);

  /* loop over all local trees */
  prev_offset = 0;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
#ifdef P4EST_ENABLE_DEBUG
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
#endif
    removed = 0;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into coarsen tree %lld with %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);

    /* state information */
    window = 0;                 /* start position of sliding window in array */
    start = 1;                  /* start position of hole in window/array */
    length = 0;                 /* length of hole in window/array */

    /* run through the array and coarsen recursively */
    incount = tquadrants->elem_count;
    while (window + P4EST_CHILDREN + length <= incount) {
      P4EST_ASSERT (window < start);

      cidz = incount;
      isfamily = 1;
      for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
        c[zz] = (window + zz < start) ?
          p4est_quadrant_array_index (tquadrants, window + zz) :
          p4est_quadrant_array_index (tquadrants, window + length + zz);

        if (zz != (size_t) p4est_quadrant_child_id (c[zz])) {
          isfamily = 0;
          if (callback_orphans) {
            c[1] = NULL;
            (void) coarsen_fn (p4est, jt, c);
          }
          break;
        }
      }
      /* in a complete tree, the only way P4EST_CHILDREN consecutive quadrants
       * have the correct consecutive child_id's is if they are, in fact, a
       * family.
       */
      P4EST_ASSERT (!isfamily || p4est_quadrant_is_familypv (c));
      if (isfamily && coarsen_fn (p4est, jt, c)) {
        /* coarsen this family of quadrants */
        if (replace_fn == NULL) {
          for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
            p4est_quadrant_free_data (p4est, c[zz]);
          }
        }
        tree->quadrants_per_level[c[0]->level] -= P4EST_CHILDREN;
        cfirst = c[0];
        if (replace_fn != NULL) {
          qtemp = *(c[0]);
          c[0] = &qtemp;
        }
        p4est_quadrant_parent (c[0], cfirst);
        p4est_quadrant_init_data (p4est, jt, cfirst, init_fn);
        tree->quadrants_per_level[cfirst->level] += 1;
        p4est->local_num_quadrants -= P4EST_CHILDREN - 1;
        removed += P4EST_CHILDREN - 1;

        cidz = (size_t) p4est_quadrant_child_id (cfirst);
        start = window + 1;
        length += P4EST_CHILDREN - 1;

        if (replace_fn != NULL) {
          replace_fn (p4est, jt, P4EST_CHILDREN, c, 1, &cfirst);
          for (zz = 0; zz < P4EST_CHILDREN; zz++) {
            p4est_quadrant_free_data (p4est, c[zz]);
          }
        }
      }

      if (cidz <= window && coarsen_recursive) {
        window -= cidz;
      }
      else {
        ++window;
        if (window == start && start + length < incount) {
          if (length > 0) {
            cfirst = p4est_quadrant_array_index (tquadrants, start);
            clast = p4est_quadrant_array_index (tquadrants, start + length);
            *cfirst = *clast;
          }
          start = window + 1;
        }
      }
    }

    /* adjust final array size */
    if (length > 0) {
      for (zz = start + length; zz < incount; ++zz) {
        cfirst = p4est_quadrant_array_index (tquadrants, zz - length);
        clast = p4est_quadrant_array_index (tquadrants, zz);
        *cfirst = *clast;
      }
      sc_array_resize (tquadrants, incount - length);
    }

    /* call remaining orphans */
    if (callback_orphans) {
      c[1] = NULL;
      for (zz = window; zz < incount - length; ++zz) {
        c[0] = p4est_quadrant_array_index (tquadrants, zz);
        (void) coarsen_fn (p4est, jt, c);
      }
    }

    /* compute maximum level */
    maxlevel = 0;
    num_quadrants = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
      num_quadrants += tree->quadrants_per_level[i];    /* same type */
      if (tree->quadrants_per_level[i] > 0) {
        maxlevel = i;
      }
    }
    tree->maxlevel = (int8_t) maxlevel;
    tree->quadrants_offset = prev_offset;
    prev_offset += num_quadrants;

    /* do some sanity checks */
    P4EST_ASSERT (num_quadrants == (p4est_locidx_t) tquadrants->elem_count);
    P4EST_ASSERT (tquadrants->elem_count == incount - removed);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size - removed ==
                    p4est->user_data_pool->elem_count);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done coarsen tree %lld now %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);
  }
  if (p4est->last_local_tree >= 0) {
    for (; jt < p4est->connectivity->num_trees; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  P4EST_FUNC_SHOT (p4est, &snap);
}

/** Check if the insulation layer of a quadrant overlaps anybody.
 * If yes, the quadrant itself is scheduled for sending.
 * Both quadrants are in the receiving tree's coordinates.
 * \param [in]  qtree       Tree id of the receiving tree.
 * \param [in]  inter_tree  Boolean flag to specify inter-tree communication.
 * \param [in]  q           The quadrant to be sent if there is overlap.
 * \param [in]  insul       An insulation quadrant of \a q.
 * \param [in,out]  first_peer  Lowest peer, will be updated.
 * \param [in,out]  last_peer   Highest peer, will be updated.
 */
static void
p4est_balance_schedule (p4est_t * p4est, p4est_balance_peer_t * peers,
                        p4est_topidx_t qtree, int inter_tree,
                        const p4est_quadrant_t * q,
                        const p4est_quadrant_t * insul,
                        int *first_peer, int *last_peer)
{
  const int           rank = p4est->mpirank;
  int                 found;
  int                 back, pos;
  int                 owner, first_owner, last_owner;
  p4est_gloidx_t     *global_first_quadrant = p4est->global_first_quadrant;
  p4est_quadrant_t    ld, *s;
  p4est_balance_peer_t *peer;

  P4EST_QUADRANT_INIT (&ld);

  /* querying insul is equivalent to querying first descendant */
  first_owner = p4est_comm_find_owner (p4est, qtree, insul, rank);
  /* querying last descendant */
  p4est_quadrant_last_descendant (insul, &ld, P4EST_QMAXLEVEL);
  last_owner = p4est_comm_find_owner (p4est, qtree, &ld, rank);

  /* send to all processors possibly intersecting insulation */
  for (owner = first_owner; owner <= last_owner; ++owner) {
    if (owner == rank && !inter_tree) {
      /* do not send to self for the same tree */
      continue;
    }
    if (global_first_quadrant[owner] == global_first_quadrant[owner + 1]) {
      /* do not send to empty processors */
      continue;
    }
    peer = peers + owner;
    /* avoid duplicates in the send array */
    found = 0;
    for (back = 0; back < P4EST_INSUL - 1; ++back) {
      pos = (int) peer->send_first.elem_count - back - 1;
      if (pos < 0) {
        break;
      }
      s = (p4est_quadrant_t *) sc_array_index_int (&peer->send_first, pos);
      if (p4est_quadrant_is_equal (s, q) && s->p.piggy2.which_tree == qtree
          && s->p.piggy2.from_tree == q->p.piggy2.from_tree
          && s->pad16 == q->pad16) {
        found = 1;
        break;
      }
    }
    if (found) {
      continue;
    }

    /* copy quadrant into shipping list */
    s = p4est_quadrant_array_push (&peer->send_first);
    *s = *q;
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    s->p.piggy2.which_tree = qtree;     /* piggy back tree id */

    /* update lowest and highest peer */
    if (owner != rank) {
      *first_peer = SC_MIN (owner, *first_peer);
      *last_peer = SC_MAX (owner, *last_peer);
    }
  }
}

static void
p4est_balance_response (p4est_t * p4est, p4est_balance_peer_t * peer,
                        p4est_connect_type_t balance, sc_array_t * borders)
{
  sc_array_t         *first_seeds = sc_array_new (sizeof (p4est_quadrant_t));

/* compute and uniqify overlap quadrants */
  p4est_tree_compute_overlap (p4est, &peer->recv_first,
                              &peer->send_second, balance, borders,
                              first_seeds);
  /* replace peer->recv_first with first_seeds */
  p4est_tree_uniqify_overlap (&peer->send_second);
  p4est_tree_uniqify_overlap (first_seeds);
  /* replace peer->recv_first with first_seeds */
  sc_array_resize (&peer->recv_first, first_seeds->elem_count);
  memcpy (peer->recv_first.array, first_seeds->array,
          first_seeds->elem_size * first_seeds->elem_count);
  sc_array_destroy (first_seeds);

  if (p4est->inspect) {
    p4est->inspect->balance_comm_sent += peer->send_second.elem_count;
    if (peer->send_second.elem_count) {
      p4est->inspect->balance_comm_nzpeers++;
    }
  }
}

void
p4est_balance (p4est_t * p4est, p4est_connect_type_t btype,
               p4est_init_t init_fn)
{
  p4est_balance_ext (p4est, btype, init_fn, NULL);
}

typedef struct
{
  p4est_topidx_t nt;
  int nl;
  int u[P4EST_UTRANSFORM];
}
p4est_tree_neigh_t;

typedef struct
{
  int offset[P4EST_INSUL + 1];
  sc_array_t tnarray;
}
p4est_tree_neigh_info_t;

const static int p4est_insul_lookup[P4EST_INSUL] = {
#ifndef P4_TO_P8
  0,  2,  1,  0, -1,  1,  2,  3,  3,
#else
  0,  0,  1,  4,  4,  5,  2,  1,  3,
  8,  2,  9,  0, -1,  1, 10,  3, 11,
  4,  2,  5,  6,  5,  7,  6,  3,  7,
#endif
};
const static int p4est_face_to_insul[P4EST_FACES] = {
#ifndef P4_TO_P8
  3,  5,  1,  7,
#else
  12, 14, 10, 16, 4, 22,
#endif
};
#ifdef P4_TO_P8
const static int p8est_edge_to_insul[P8EST_EDGES] = { 1, 7, 19, 25, 3, 5, 21, 23, 9, 11, 15, 17 };
#endif
const static int p4est_corner_to_insul[P4EST_CHILDREN] = {
#ifndef P4_TO_P8
  0, 2, 6, 8,
#else
  0, 2, 6, 8, 18, 20, 24, 26,
#endif
};

static int
p4est_tree_get_conn_info (p4est_t *p4est, p4est_topidx_t t,
                                  p4est_tree_neigh_info_t *info)
{
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t nt;
  p4est_quadrant_t root, neigh;
  int i, j, k, l, nl;
  int kstart, kend;
  size_t s, nneigh;
  p4est_corner_info_t ci;
  sc_array_t *cta = &ci.corner_transforms;
#ifdef P4_TO_P8
  p8est_edge_info_t ei;
  sc_array_t *eta = &ei.edge_transforms;
#endif

#ifdef P4_TO_P8
  kstart = 0;
  kend = 3;
#else
  kstart = 1;
  kend = 2;
#endif

  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#ifdef P4_TO_P8
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif

  sc_array_init_size (&(info->tnarray), sizeof (p4est_tree_neigh_t), P4EST_INSUL);
  sc_array_truncate (&(info->tnarray));

  memset (&root, 0, sizeof (p4est_quadrant_t));

  /* loop over the tree insulation layer,
   * cache neighboring trees and transforms */
  for (k = kstart, l = 0; k < kend; k++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++, l++) {
        int type = (k != 1) + (j != 1) + (i != 1);
        int id = p4est_insul_lookup[l];

        P4EST_ASSERT (l < P4EST_INSUL);

        info->offset[l] = (int) info->tnarray.elem_count;

        if (!type) {
          p4est_tree_neigh_t *tn = sc_array_push (&info->tnarray);

          tn->nt = t;
          tn->nl = P4EST_INSUL / 2;
        }
        else if (type == 1) {
          /* outside face */
          nt = p4est_quadrant_face_neighbor_extra (&root, t, id, &neigh, &nl, conn);
          nl = nl % P4EST_FACES;
          if (nt >= 0) {
            int f[P4EST_FTRANSFORM];
            p4est_tree_neigh_t *tn = (p4est_tree_neigh_t *) sc_array_push (&(info->tnarray));

            tn->nt = nt;
            tn->nl = p4est_face_to_insul[nl];
            p4est_find_face_transform (conn, t, id, f);
            p4est_face_transform_to_utransform (f, &(tn->u[0]));
          }
        }
#ifdef P4_TO_P8
        else if (type == 2) {
          /* outside edge */
          p4est_tree_neigh_t *tn;

          p8est_find_edge_transform (conn, t, id, &ei);
          nneigh = eta->elem_count;
          tn = (p4est_tree_neigh_t *) sc_array_push_count (&(info->tnarray), nneigh);
          for (s = 0; s < nneigh; s++) {
            p8est_edge_transform_t *et = p8est_edge_array_index (eta, s);

            tn[s].nt = et->ntree;
            tn[s].nl = p8est_edge_to_insul[et->nedge];
            p8est_edge_transform_to_utransform (et, id, &(tn[s].u[0]));
          }
          sc_array_truncate (eta);
        }
#endif
        else {
          /* outside corner */
          p4est_tree_neigh_t *tn;

          p4est_find_corner_transform (conn, t, id, &ci);
          nneigh = cta->elem_count;
          tn = (p4est_tree_neigh_t *) sc_array_push_count (&(info->tnarray), nneigh);
          for (s = 0; s < nneigh; s++) {
            p4est_corner_transform_t *ct = p4est_corner_array_index (cta, s);

            tn[s].nt = ct->ntree;
            tn[s].nl = p4est_corner_to_insul[ct->ncorner];
            p4est_corner_transform_to_utransform (ct, id, &(tn[s].u[0]));
          }
        }
      }
    }
  }
  info->offset[P4EST_INSUL] = info->tnarray.elem_count;

  sc_array_reset (cta);
#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  return 0;
}

static int
p4est_root_insul (p4est_quadrant_t *q)
{
  p4est_qcoord_t x[P4EST_DIM];
  int i;
  int total;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  x[0] = (q->x >> P4EST_MAXLEVEL);
  x[1] = (q->y >> P4EST_MAXLEVEL);
#ifdef P4_TO_P8
  x[2] = (q->z >> P4EST_MAXLEVEL);
#endif
  for (i = P4EST_DIM - 1, total = 0; i >= 0; i--) {
    total *= 3;
    P4EST_ASSERT (-1 <= x[i] && x[i] <= 1);
    total += x[i] + 1;
  }
  P4EST_ASSERT (0 <= total && total < P4EST_INSUL);
  return total;
}

static void
p4est_tree_ghost_neigh_add_r (p4est_t *p4est, p4est_topidx_t from_tree, int insul, int insul_max, p4est_topidx_t which_tree, int nlevel,
                              int32_t touch, int p_fi, int p_la, p4est_quadrant_t *fi, p4est_quadrant_t *la,
                              int *proc_hash, sc_array_t *procs)
{
  int32_t rb;
  p4est_quadrant_t fi_stack[64];
  p4est_quadrant_t la_stack[64];
  int p_fi_stack[64];
  int p_la_stack[64];
  int stack_pointer = 0;

  fi_stack[0] = *fi;
  la_stack[0] = *la;
  p_fi_stack[0] = p_fi;
  p_la_stack[0] = p_la;

  for (;;) {
    P4EST_ASSERT (stack_pointer < 64);
    rb = p4est_find_range_boundaries (&fi_stack[stack_pointer], &la_stack[stack_pointer], nlevel, NULL,
#ifdef P4_TO_P8
                                      NULL,
#endif
                                      NULL);
    if (!(rb & touch)) {
      /* move down the stack */
      while (stack_pointer > 0 && p_la_stack[stack_pointer] == p_la_stack[stack_pointer - 1]) {
        stack_pointer--;
      }
      if (!stack_pointer) {
        return;
      }
      p_fi_stack[stack_pointer] = p_fi_stack[stack_pointer - 1];
      p_la_stack[stack_pointer] = p_la_stack[stack_pointer - 1];
      fi_stack[stack_pointer] = fi_stack[stack_pointer - 1];
      la_stack[stack_pointer] = la_stack[stack_pointer - 1];
      continue;
    }
    if (p_fi_stack[stack_pointer] == p_la_stack[stack_pointer]) { /* (uniquely) add this tree's portion of the process range to the neighborhood */
      int p = p_fi_stack[stack_pointer];
      int idx = proc_hash[p];

      if (!idx) {
        /* process has not been seen before, add it */
        *((int *) sc_array_push (procs)) = p;
        proc_hash[p] = (int) procs->elem_count;
      }
      while (stack_pointer > 0 && p_la_stack[stack_pointer] == p_la_stack[stack_pointer - 1]) {
        stack_pointer--;
      }
      if (!stack_pointer) {
        return;
      }
      p_fi_stack[stack_pointer] = p_fi_stack[stack_pointer - 1];
      p_la_stack[stack_pointer] = p_la_stack[stack_pointer - 1];
      fi_stack[stack_pointer] = fi_stack[stack_pointer - 1];
      la_stack[stack_pointer] = la_stack[stack_pointer - 1];
      continue;
    }
    /* recurse only if this range of processes is adjacent */
    {
      int p = p_fi_stack[stack_pointer];
      int q = p_la_stack[stack_pointer];
      int num_procs = q + 1 - p;
      int mid = p + num_procs / 2;
      uint64_t id;
      p4est_quadrant_t *lalo;

      P4EST_ASSERT (stack_pointer < 63);
      p_fi_stack[stack_pointer + 1] = p;
      p_la_stack[stack_pointer + 1] = mid - 1;
      fi_stack[stack_pointer + 1] = fi_stack[stack_pointer];

      lalo = &p4est->global_first_position[mid];
      id = p4est_quadrant_linear_id (lalo, P4EST_QMAXLEVEL);
      p4est_quadrant_set_morton (&la_stack[stack_pointer + 1], P4EST_QMAXLEVEL, id - 1);
      p_fi_stack[stack_pointer] = mid;
      fi_stack[stack_pointer] = *lalo;

      stack_pointer++;
    }
  }
}

static void
p4est_tree_ghost_neighborhood_add (p4est_t *p4est, int from_tree, int insul, int insul_max, p4est_topidx_t which_tree, p4est_quadrant_t *q,
                                   p4est_quadrant_t *n, int *proc_hash, sc_array_t *procs)
{
  p4est_qcoord_t diff[3];
  int ninsul;
  int hshift;
  int i;
  int l = q->level;
  int32_t touch_table[P4EST_INSUL] = {
#ifndef P4_TO_P8
    0x10, 0x04, 0x20,
    0x01,   -1, 0x02,
    0x40, 0x08, 0x80,
#else
    0x0040000, 0x0000040, 0x0080000,
    0x0000400, 0x0000010, 0x0000800,
    0x0100000, 0x0000080, 0x0200000,

    0x0004000, 0x0000004, 0x0008000,
    0x0000001,        -1, 0x0000002,
    0x0010000, 0x0000008, 0x0020000,

    0x0400000, 0x0000100, 0x0800000,
    0x0001000, 0x0000020, 0x0002000,
    0x1000000, 0x0000200, 0x2000000,
#endif
  };
  int32_t touch;
  p4est_quadrant_t fi, la;
  int rank = p4est->mpirank;
  int p_fi, p_la;

  P4EST_ASSERT (n->level == l);
  P4EST_ASSERT (p4est_quadrant_is_valid (n));
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  hshift = P4EST_MAXLEVEL - l;
  diff[0] = (q->x - n->x) >> hshift;
  diff[1] = (q->y - n->y) >> hshift;
#ifdef P4_TO_P8
  diff[2] = (q->z - n->z) >> hshift;
#else
  diff[2] = 0;
#endif
  ninsul = 0;
  for (i = P4EST_DIM - 1, ninsul = 0; i >= 0; i--) {
    ninsul *= 3;
    P4EST_ASSERT (-1 <= diff[i] && diff[i] <= 1);
    ninsul += diff[i] + 1;
  }
  P4EST_ASSERT (ninsul != P4EST_INSUL / 2);
  touch = touch_table[ninsul];

  p4est_quadrant_first_descendant (n, &fi, P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (n, &la, P4EST_QMAXLEVEL);
  p_fi = p4est_comm_find_owner (p4est, which_tree, &fi, rank);
  p_la = p4est_comm_find_owner (p4est, which_tree, &la, p_fi);
  if (p_fi == rank && p_la == rank) {
    return;
  }
  p4est_tree_ghost_neigh_add_r (p4est, from_tree, insul, insul_max, which_tree, n->level, touch, p_fi, p_la, &fi, &la, proc_hash, procs);
}

static void
p4est_tree_ghost_neighborhood_insert (p4est_t *p4est, p4est_topidx_t t, p4est_quadrant_t *q,
                                      p4est_tree_neigh_info_t *info, int *proc_hash, sc_array_t *procs)
{
  p4est_qcoord_t h = P4EST_QUADRANT_LEN (q->level);
  int i, j, k, l, insulmax;
#ifdef P4_TO_P8
  int kstart = 0, kend = 3;
#else
  int kstart = 1, kend = 2;
#endif

  insulmax = info->offset[P4EST_INSUL];
  for (k = kstart, l = 0; k < kend; k++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++, l++) {
        p4est_quadrant_t n;
        p4est_quadrant_t nn, nq;
        int rootinsul;
        size_t s;

        if (l == P4EST_INSUL / 2) {
          continue;
        }

        n = *q;
        n.x += (i - 1) * h;
        n.y += (j - 1) * h;
#ifdef P4_TO_P8
        n.z += (k - 1) * h;
#endif

        rootinsul = p4est_root_insul (&n);
        if (rootinsul == P4EST_INSUL / 2) {
          p4est_tree_ghost_neighborhood_add (p4est, t, info->offset[rootinsul], insulmax, t, q, &n, proc_hash, procs);
        }
        else {
          size_t offset = info->offset[rootinsul];
          size_t stop = info->offset[rootinsul + 1];

          for (s = offset; s < stop; s++) {
            p4est_tree_neigh_t *tn = (p4est_tree_neigh_t *) sc_array_index (&info->tnarray, s);

            p4est_quadrant_utransform (q, &nq, &(tn->u[0]), 0);
            p4est_quadrant_utransform (&n, &nn, &(tn->u[0]), 0);
            p4est_tree_ghost_neighborhood_add (p4est, t, s, insulmax, tn->nt, &nq, &nn, proc_hash, procs);
          }
        }
      }
    }
  }
}

static void
p4est_tree_ghost_neighborhood (p4est_t *p4est, p4est_topidx_t t, int *proc_hash, sc_array_t *procs, p4est_tree_neigh_info_t *info)
{
  p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, t);
  p4est_quadrant_t f = tree->first_desc;
  p4est_quadrant_t l = tree->last_desc;
  p4est_quadrant_t a, af, al, fstop, temp;
  int           alevel, fid, lid, idxor, dir, fstopid;

  p4est_nearest_common_ancestor (&f, &l, &a);

  alevel = a.level;

  p4est_quadrant_first_descendant (&a, &af, P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (&a, &al, P4EST_QMAXLEVEL);
  if (p4est_quadrant_is_equal (&f, &af) &&
      p4est_quadrant_is_equal (&l, &al)) {
    p4est_tree_ghost_neighborhood_insert (p4est, t, &a, info, proc_hash, procs);
    return;
  }
  fid = p4est_quadrant_ancestor_id (&f, alevel + 1);
  lid = p4est_quadrant_ancestor_id (&l, alevel + 1);
  P4EST_ASSERT (lid > fid);
  idxor = lid ^ fid;
  dir = SC_LOG2_8 (idxor);
  P4EST_ASSERT (dir >= 0);
  fstopid = (lid >> dir) << dir;
  p4est_quadrant_child (&a, &fstop, fstopid);
  while (f.level > alevel + 1 &&
         p4est_quadrant_child_id (&f) == 0) {
    p4est_quadrant_parent (&f, &f);
  }
  while (l.level > alevel + 1 &&
         p4est_quadrant_child_id (&l) == P4EST_CHILDREN - 1) {
    p4est_quadrant_parent (&l, &l);
  }
  P4EST_ASSERT (!p4est_quadrant_overlaps (&f, &fstop) && p4est_quadrant_compare (&f, &fstop) < 0);
  P4EST_ASSERT (p4est_quadrant_compare (&fstop, &l) <= 0);
  for (;;) {
    p4est_tree_ghost_neighborhood_insert (p4est, t, &f, info, proc_hash, procs);
    fid = p4est_quadrant_child_id (&f);
    while (fid == P4EST_CHILDREN - 1) {
      p4est_quadrant_parent (&f, &f);
      fid = p4est_quadrant_child_id (&f);
    }
    p4est_quadrant_sibling (&f, &temp, fid + 1);
    f = temp;
    if (p4est_quadrant_is_equal (&f, &fstop)) {
      break;
    }
  }
  for (;;) {
    while (p4est_quadrant_is_ancestor (&fstop, &l)) {
      fstop.level++;
    }
    p4est_tree_ghost_neighborhood_insert (p4est, t, &fstop, info, proc_hash, procs);
    if (p4est_quadrant_is_equal (&fstop, &l)) {
      break;
    }
    fstopid = p4est_quadrant_child_id (&fstop);
    P4EST_ASSERT (fstopid < P4EST_CHILDREN - 1);
    p4est_quadrant_sibling (&fstop, &temp, fstopid + 1);
    fstop = temp;
  }
}

static int
p4est_balance_sort_min_insulation_level (p4est_t *p4est, p4est_tree_neigh_info_t *info, p4est_quadrant_t *q)
{
  int i, j, k, l, m, p;
#ifdef P4_TO_P8
  int kstart = 0, kend = 3;
#else
  int kstart = 1, kend = 2;
#endif
  p4est_quadrant_t fd, ld;
  p4est_quadrant_t *nd;

  for (l = 0; l < P4EST_QMAXLEVEL; l++) {
    p4est_quadrant_t a, n;
    p4est_qcoord_t h;

    p4est_quadrant_ancestor (q, l, &a);
    h = P4EST_QUADRANT_LEN (l);
    n.level = l;
    n.p.which_tree = q->p.which_tree;
    for (m = 0, k = kstart; k < kend; k++) {
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++, m++) {
          int rootinsul;

          n.x = a.x + (i - 1) * h;
          n.y = a.y + (j - 1) * h;
#ifdef P4_TO_P8
          n.z = a.z + (j - 1) * h;
#endif

          rootinsul = p4est_root_insul (&n);

          if (rootinsul == P4EST_INSUL / 2) {
            /* in root */
            p4est_quadrant_first_descendant (&n, &fd, P4EST_QMAXLEVEL);
            p4est_quadrant_last_descendant (&n, &ld, P4EST_QMAXLEVEL);
            ld.p.which_tree = n.p.which_tree;
            p = p4est_comm_find_owner (p4est, q->p.which_tree, &fd, p4est->mpirank);
            nd = &p4est->global_first_position[p + 1];
            if (p4est_quadrant_compare_piggy (&ld, nd) < 0) {
              /* n is entirely contained in one process */
              return l;
            }
          }
          else {
            size_t s, s_start, s_end;

            s_start = info->offset[rootinsul];
            s_end = info->offset[rootinsul + 1];
            for (s = s_start; s < s_end; s++) {
              p4est_quadrant_t nn;
              p4est_tree_neigh_t *tn = (p4est_tree_neigh_t *) sc_array_index (&(info->tnarray), s);

              p4est_quadrant_utransform (&n, &nn, &(tn->u[0]), 0);
              p4est_quadrant_first_descendant (&nn, &fd, P4EST_QMAXLEVEL);
              p4est_quadrant_last_descendant (&nn, &ld, P4EST_QMAXLEVEL);
              ld.p.which_tree = tn->nt;
              p = p4est_comm_find_owner (p4est, tn->nt, &fd, p4est->mpirank);
              nd = &p4est->global_first_position[p + 1];
              if (p4est_quadrant_compare_piggy (&ld, nd) < 0) {
                /* nn is entirely contained in one process */
                return l;
              }
            }
          }
        }
      }
    }
  }
  SC_ABORT_NOT_REACHED();
  return -1;
}

static void
sc_array_crop (sc_array_t *array, size_t a_start, size_t a_end)
{
  P4EST_ASSERT (array);
  P4EST_ASSERT (a_end >= a_start);
  P4EST_ASSERT (a_end <= array->elem_count);
  if (a_end == a_start) {
    array->elem_count = 0;
    return;
  }
  if (a_start) {
    memmove (array->array, sc_array_index (array, a_start), (a_end - a_start) * array->elem_size);
  }
  array->elem_count = (a_end - a_start);
}

static void
sc_array_push_copy (sc_array_t *array, size_t a_start, size_t a_end, sc_array_t *result)
{
  P4EST_ASSERT (result && array != result);
  P4EST_ASSERT (result->elem_size == array->elem_size);
  P4EST_ASSERT (a_end >= a_start);

  if (a_end > a_start) {
    char *dest = sc_array_push_count (result, a_end - a_start);
    memcpy (dest, sc_array_index (array, a_start), (a_end - a_start) * array->elem_size);
  }
}

static void
sc_array_push_cut (sc_array_t *array, size_t a_start, size_t a_end, sc_array_t *result)
{
  P4EST_ASSERT (result && array != result);
  P4EST_ASSERT (result->elem_size == array->elem_size);
  P4EST_ASSERT (a_end >= a_start);

  if (a_end > a_start) {
    char *dest = sc_array_push_count (result, a_end - a_start);

    memcpy (dest, sc_array_index (array, a_start), (a_end - a_start) * array->elem_size);
    if (a_end < array->elem_count) {
      memmove (sc_array_index (array, a_start), sc_array_index (array, a_end), (array->elem_count - a_end) * array->elem_size);
    }
    array->elem_count -= (a_end - a_start);
  }
}

#if 0
static void
p4est_quadrant_compare_insertion (const void *key, const void *memb)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) key;
  const p4est_quadrant_t *p = (const p4est_quadrant_t *) memb;

  if (p4est_quadrant_compare (q, &p[0]) < 0) {
    return -1;
  }
  if (p4est_quadrant_compare (q, &p[1]) < 0) {
    return 0;
  }
  return 1;
}

static size_t
p4est_find_insertion (sc_array_t *qarray, p4est_quadrant_t *key)
{
  p4est_quadrant_t *q;
  P4EST_ASSERT (qarray->elem_size == sizeof (p4est_quadrant_t));

  if (!(qarray->elem_count)) {
    return 0;
  }
  q = (p4est_quadrant_t *) qarray->array;
  if (p4est_quadrant_compare (key, q) < 0) {
    return 0;
  }
  else {
    size_t sz;
    sc_array_t comparator;

    sc_array_init_view (&comparator, qarray, 0, qarray->elem_count - 1);
    sz = sc_array_bsearch (&comparator, key, p4est_quadrant_compare_insertion);
    if (sz < 0) {
      return qarray->elem_count;
    }
    else {
      return (size_t) sz;
    }
  }
}
#endif

enum {
  P4EST_BALSORT_ISEND,
  P4EST_BALSORT_IPROBE,
  P4EST_BALSORT_IRECV,
  P4EST_BALSORT_MERGE,
  P4EST_BALSORT_WAIT,
  P4EST_BALSORT_SYNC,
};

static void
p4est_quadrant_array_merge_reduce (sc_array_t *A, sc_array_t *B)
{
  size_t nnew = B->elem_count;
  size_t nold = A->elem_count;
  p4est_quadrant_t *a, *b, *c;
  p4est_quadrant_t cp;
  size_t i, j, k;

  if (!nnew) {
    return;
  }
  P4EST_ASSERT (nold == 0 || p4est_quadrant_array_is_reduced (A));
  P4EST_ASSERT (p4est_quadrant_array_is_reduced (B));
  (void) sc_array_push_count (A, nnew);
  if (!nold) {
    sc_array_copy (A, B);
    return;
  }
  a = (p4est_quadrant_t *) sc_array_index (A, nnew);
  b = (p4est_quadrant_t *) B->array;
  c = (p4est_quadrant_t *) A->array;
  memmove (a, c, nold * A->elem_size);
  for (i = 0, j = 0, k = 0; i + j < nold + nnew;) {
    p4est_quadrant_t *q, qp;
    int comp;

    if (i == nold) {
      comp = -1;
    }
    else if (j == nnew) {
      comp = 1;
    }
    else {
      comp = p4est_quadrant_compare (b, a);
    }
    if (comp < 0) {
      q = b++;
      j++;
    }
    else if (comp > 0) {
      q = a++;
      i++;
    }
    else {
      q = a++;
      b++;
      i++;
      j++;
    }
    if (!q->level) {
      continue;
    }
    if (!k) {
      p4est_quadrant_parent (q, &cp);
      *c = *q;
      k++;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_compare (c, q) < 0);
      p4est_quadrant_parent (q, &qp);
      if (p4est_quadrant_overlaps (&cp, &qp)) {
        if (q->level > c->level) {
          *c = *q;
          cp = qp;
        }
      }
      else {
        *(++c) = *q;
        cp = qp;
        k++;
      }
    }
  }
  A->elem_count = k;
  P4EST_ASSERT (p4est_quadrant_array_is_reduced (A));
}

/* divide up a sorted list of quadrants from a given tree between the
 * neighbor processes */
static void
p4est_balance_sort_divide (p4est_t *p4est, sc_array_t *procs, sc_array_t *proc_bufs,
                           p4est_topidx_t which_tree, sc_array_t *tquads, int skip_self)
{
  size_t n, nneigh = procs->elem_count;

  for (n = 0; n < nneigh; n++) {
    int p = *((int *) sc_array_index (procs, n));
    p4est_quadrant_t fp, np, lb, ub;
    p4est_quadrant_t *nq;
    ssize_t sz;
    size_t qstart, qend;
    size_t nbuf, s;

    if (p == p4est->mpirank && skip_self) {
      continue;
    }

    fp = p4est->global_first_position[p];
    np = p4est->global_first_position[p+1];
    if (np.p.which_tree < which_tree || fp.p.which_tree > which_tree) {
      continue;
    }
    if (fp.p.which_tree == which_tree) {
      lb = fp;

      while (lb.level && !p4est_quadrant_child_id(&lb)) {
        lb.level--;
      }
      sz = p4est_find_lower_bound (tquads, &lb, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qstart = sz;
    }
    else {
      qstart = 0;
    }

    if (np.p.which_tree > which_tree) {
      qend = tquads->elem_count;
    }
    else {
      ub = np;

      while (ub.level && !p4est_quadrant_child_id (&ub)) {
        ub.level--;
      }
      sz = p4est_find_lower_bound (tquads, &ub, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qend = sz;
    }
    P4EST_ASSERT (qend >= qstart);
    if (qend > qstart) {
      nbuf = proc_bufs[n].elem_count;
      sc_array_push_copy (tquads, qstart, qend, &proc_bufs[n]);
      nq = p4est_quadrant_array_index (&proc_bufs[n], nbuf);
      for (s = 0; s < qend - qstart; s++) {
        nq[s].pad8 = 0;
        nq[s].pad16 = 0;
        nq[s].p.user_data = NULL;
        nq[s].p.which_tree = which_tree;
      }
    }
  }
}

static void
p4est_balance_sort_compute_pattern (p4est_t *p4est,
                                    int flt, int llt, int num_trees,
                                    p4est_tree_neigh_info_t *tinfo,
                                    int *proc_hash, sc_array_t *procs,
                                    int (*procrange)[2], int *minlevel, int *num_sort,
                                    p4est_quadrant_t *desc)
{
  p4est_topidx_t t;
  sc_flopinfo_t snap;
  size_t s;
  int rank = p4est->mpirank;

  P4EST_FUNC_SNAP (p4est, &snap);

  /* first figure out the parallel neighborhood */
  sc_array_truncate (procs);
  *((int *) sc_array_push (procs)) = rank;
  proc_hash[rank] = 1;
  /* block balance local */
  for (t = flt; t <= llt; t++) {
    p4est_tree_neigh_info_t *info;

    /* Get the connectivity information for this tree */
    info = &tinfo[t - flt];
    p4est_tree_get_conn_info (p4est, t, info);

    /* find ghost neighbors of this tree */
    if (num_trees) {
      p4est_tree_ghost_neighborhood (p4est, t, proc_hash, procs, info);
    }
  }
  sc_array_sort (procs, sc_int_compare);
  for (s = 0; s < procs->elem_count; s++) {
    int p = *((int *) sc_array_index (procs, s));

    proc_hash[p] = s + 1;
  }
#ifdef P4EST_ENABLE_DEBUG
  {
    int *proc_hash_T = P4EST_ALLOC (int, p4est->mpisize);
    int mpiret;
    int i;

    mpiret = sc_MPI_Alltoall (proc_hash, 1, sc_MPI_INT, proc_hash_T, 1, sc_MPI_INT,
                              p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < p4est->mpisize; i++) {
      P4EST_ASSERT (!!proc_hash[i] == !!proc_hash_T[i]);
    }
    P4EST_FREE (proc_hash_T);
  }
#endif

  /* compute the processes I communicate with in the sort phase */
  if (llt >= flt) {
    int i;
    int empty;

    desc[0] = p4est->global_first_position[rank];
    desc[1] = p4est->global_first_position[rank + 1];
    empty = p4est_quadrant_is_equal_piggy (&desc[0], &desc[1]);
    if (desc[1].p.which_tree > llt) {
      p4est_qcoord_t last = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
      desc[1].p.which_tree--;

      desc[1].x = last;
      desc[1].y = last;
#ifdef P4_TO_P8
      desc[1].z = last;
#endif
    }
    else {
      uint64_t id = p4est_quadrant_linear_id (&desc[1], P4EST_QMAXLEVEL);

      p4est_quadrant_set_morton (&desc[1], P4EST_QMAXLEVEL, id - 1);
    }

    for (i = 0; i < 2; i++) {
      p4est_quadrant_t a, fd, ld;
      int t = i ? llt : flt;

      minlevel[i] = p4est_balance_sort_min_insulation_level (p4est, &tinfo[t - flt], &desc[i]);
      p4est_quadrant_ancestor (&desc[i], minlevel[i], &a);
      p4est_quadrant_first_descendant(&a, &fd, P4EST_QMAXLEVEL);
      p4est_quadrant_last_descendant(&a, &ld, P4EST_QMAXLEVEL);
      procrange[i][0] = p4est_comm_find_owner (p4est, t, &fd, rank);
      procrange[i][1] = p4est_comm_find_owner (p4est, t, &ld, rank);
    }
    if (empty) {
      *num_sort = 1;
      if (procrange[0][0] != procrange[1][0]) {
        /* communicate with no one */
        procrange[0][0] = procrange[1][0] = -1;
        procrange[0][1] = procrange[1][1] = -2;
        *num_sort = 0;
      }
    }
    else {
      if (procrange[0][0] == procrange[1][0] &&
          procrange[0][1] == procrange[1][1]) {
        *num_sort = 1;
      }
      else {
        *num_sort = 2;
      }
      P4EST_ASSERT (procrange[0][0] <= rank);
      P4EST_ASSERT (rank <= procrange[1][1]);
      if (procrange[0][1] > procrange[1][0]) {
        P4EST_ASSERT (procrange[0][0] == procrange[1][0]);
        P4EST_ASSERT (procrange[0][1] == procrange[1][1]);
      }
      else {
        P4EST_ASSERT (procrange[0][1] == rank);
        P4EST_ASSERT (procrange[1][0] == rank);
      }
    }
  }
  else {
    /* communicate with no one */
    procrange[0][0] = procrange[1][0] = -1;
    procrange[0][1] = procrange[1][1] = -2;
    num_sort = 0;
  }
#ifdef P4EST_ENABLE_DEBUG
  {
    int * comm_out = P4EST_ALLOC_ZERO (int, p4est->mpisize);
    int * comm_in = P4EST_ALLOC (int, p4est->mpisize);
    int i, j, mpiret;

    for (i = 0; i < 2; i++) {
      for (j = procrange[i][0]; j <= procrange[i][1]; j++) {
        comm_out[j] = 1;
      }
    }
    mpiret = sc_MPI_Alltoall (comm_out, 1, sc_MPI_INT, comm_in, 1, sc_MPI_INT,
                              p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < p4est->mpisize; i++) {
      P4EST_ASSERT (comm_out[i] == comm_in[i]);
    }
    P4EST_FREE (comm_out);
    P4EST_FREE (comm_in);
  }
#endif
  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort_local (p4est_t *p4est,
                          p4est_connect_type_t btype,
                          int flt, int llt, int num_trees,
                          int num_sort, int *minlevel,
                          p4est_quadrant_t *desc,
                          sc_array_t *tree_bufs, sc_array_t *sort_bufs)
{
  p4est_topidx_t t;
  const int8_t       *pre_adapt_flags = NULL;
  p4est_locidx_t all_incount;
  int bound;
  sc_flopinfo_t snap;
  sc_mempool_t       *qpool;
  sc_mempool_t       *list_alloc;

  P4EST_FUNC_SNAP (p4est, &snap);
  if (p4est->inspect != NULL) {
    pre_adapt_flags = p4est->inspect->pre_adapt_flags;
  }

  switch (p4est_connect_type_int (btype)) {
  case 0:
    bound = 1;
    break;
  case 1:
    bound = P4EST_DIM + 1;
    break;
  case P4EST_DIM:
    bound = (1 << P4EST_DIM);
    break;
#ifdef P4_TO_P8
  case 2:
    bound = 7;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }

  qpool = p4est->quadrant_pool;
  list_alloc = sc_mempool_new (sizeof (sc_link_t));

  /* block balance local */
  for (t = flt, all_incount = 0; t < flt + num_trees; t++) {
    const int8_t       *this_pre_adapt_flags = NULL;
    p4est_tree_t       *tree;
    sc_array_t         *tquadrants;
    size_t              treecount;
    p4est_quadrant_t    root;
    size_t              ancstart, ancend, selfstart, selfend;
    sc_array_t         *inquads = &tree_bufs[t - flt];

    tree = p4est_tree_array_index (p4est->trees, t);
    tquadrants = &tree->quadrants;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into balance tree %lld with %llu\n", (long long) t,
                    (unsigned long long) tquadrants->elem_count);

    p4est_nearest_common_ancestor (&tree->first_desc, &tree->last_desc, &root);
    if (t == flt && t == llt) {
      if (minlevel[0] <= root.level) {
        P4EST_ASSERT (minlevel[1] == minlevel[0]);
        p4est_quadrant_ancestor (&tree->first_desc, minlevel[0], &root);
      }
      /* is it worth it to split the local balance when it consists of two
       * regions that do not influence each other ? */
#if 0
      else if (minlevel[0] > root.level + 1 && minlevel[1] > root.level + 1) {
        int fid, lid, level, lolevel, hilevel;
        p4est_quadrant_t lo, hi, ub, lb, temp;
        p4est_qcoord_t H, diff[P4EST_DIM];

        fid = p4est_quadrant_ancestor_id (&tree->first_desc, root.level + 1);
        lid = p4est_quadrant_ancestor_id (&tree->last_desc, root.level + 1);
        P4EST_ASSERT (lid == fid + 1);

        p4est_quadrant_child (&root, &lo, fid);
        p4est_quadrant_last_descendant (&lo, &ub, P4EST_QMAXLEVEL);
        p4est_nearest_common_ancestor (&tree->first_desc, &ub, &lo);
        lolevel = SC_MIN (lo.level, minlevel[0]);

        p4est_quadrant_child (&root, &hi, lid);
        p4est_quadrant_first_descendant (&hi, &lb, P4EST_QMAXLEVEL);
        p4est_nearest_common_ancestor (&tree->last_desc, &lb, &hi);
        hilevel = SC_MIN (hi.level, minlevel[1]);
        level = SC_MIN (lolevel, hilevel);

        p4est_quadrant_ancestor (&tree->first_desc, level, &lo);
        p4est_quadrant_ancestor (&tree->last_desc, level, &hi);

        H = P4EST_QUADRANT_LEN (level);
        diff[0] = hi.x - lo.x;
        diff[1] = hi.y - lo.y;
#ifdef P4_TO_P8
        diff[2] = hi.z - lo.z;
#endif
        if (diff[0] < -H || diff[0] > H ||
            diff[1] < -H || diff[1] > H ||
#ifdef P4_TO_P8
            diff[2] < -H || diff[2] > H ||
#endif
            0) {
          P4EST_DEBUGF ("Split root: [%d, %d], [%d, %d], [%d, %d, %d], %d\n", lolevel, hilevel, minlevel[0], minlevel[1],
                        diff[0], diff[1],
#ifdef P4_TO_P8
                        diff[2],
#else
                        0,
#endif
                        H);
        }
      }
#endif
    }
    else if (t == flt) {
      p4est_quadrant_ancestor (&tree->first_desc, SC_MIN (root.level, minlevel[0]), &root);
    }
    else if (t == llt) {
      p4est_quadrant_ancestor (&tree->first_desc, SC_MIN (root.level, minlevel[1]), &root);
    }

    if (pre_adapt_flags) {
      this_pre_adapt_flags = &pre_adapt_flags[all_incount];
    }
    all_incount += tquadrants->elem_count;


    /* local balance first pass */
    p4est_quadrant_array_reduce (tquadrants, inquads, this_pre_adapt_flags);
    if (pre_adapt_flags) {
      p4est_quadrant_array_insert_endpoints (inquads, &root, &tree->first_desc, &tree->last_desc);
    }
    p4est_balance_kernel (inquads, NULL, root.level, bound, qpool,
                          list_alloc, NULL, NULL, NULL, NULL, NULL);

    treecount = inquads->elem_count;
    P4EST_VERBOSEF ("Balance tree %lld A %llu (reduced)\n",
                    (long long) t, (unsigned long long) treecount);

    if (t > flt && t < llt) {
      continue;
    }

    if (t == flt) {
      ssize_t sz;
      p4est_quadrant_t fdesc, anc;

      fdesc = tree->first_desc;
      while (fdesc.level > 0 && !p4est_quadrant_child_id (&fdesc)) {
        fdesc.level--;
      }
      sz = p4est_find_lower_bound (inquads, &fdesc, 0);
      if (sz < 0) {
        sz = inquads->elem_count;
      }
      selfstart = sz;

      p4est_quadrant_ancestor (&desc[0], minlevel[0], &anc);

      sz = p4est_find_lower_bound (inquads, &anc, 0);
      if (sz < 0) {
        sz = inquads->elem_count;
      }
      ancstart = sz;

      P4EST_ASSERT (ancstart <= selfstart);
      sc_array_push_cut (inquads, ancstart, selfstart, &sort_bufs[0]);
      selfstart = ancstart;
    }
    else {
      selfstart = 0;
    }
    if (t == llt) {
      ssize_t sz;
      p4est_quadrant_t anc, ub;

      sz = p4est_find_higher_bound (inquads, &tree->last_desc, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      selfend = sz;

      p4est_quadrant_ancestor (&desc[1], minlevel[1], &anc);
      p4est_quadrant_last_descendant (&anc, &ub, P4EST_QMAXLEVEL);
      sz = p4est_find_higher_bound (inquads, &ub, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      ancend = sz;
      sc_array_push_cut (inquads, selfend, ancend, &sort_bufs[1]);
    }
    else {
      selfend = inquads->elem_count;
    }
#if 0
    sc_array_crop (inquads, selfstart, selfend);
#endif
  }
  sc_mempool_destroy (list_alloc);

  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort_sort (p4est_t *p4est, int flt, int llt,
                         int num_trees, int num_sort,
                         int (*procrange)[2], sc_array_t *sort_bufs)
{
  sc_flopinfo_t snap;
  int rank = p4est->mpirank;

  P4EST_FUNC_SNAP (p4est, &snap);
#ifdef P4EST_ENABLE_MPI
  if (p4est->inspect && p4est->inspect->stats) {
    if (!sc_statistics_has (p4est->inspect->stats, "balance sort send")) {
      sc_statistics_add_empty (p4est->inspect->stats, "balance sort send");
    }
    if (!sc_statistics_has (p4est->inspect->stats, "balance sort recv")) {
      sc_statistics_add_empty (p4est->inspect->stats, "balance sort recv");
    }
  }
  /* block balance parallel sort */
  if (num_sort > 0) {
    int i;
    MPI_Request sendreq[2];
    MPI_Request recvreq[2];
    int mpiret;
    int bit[2] = {0, 0};
    int proc[2] = {-1, -1};
    int state[2];
    sc_array_t recv_bufs[2], send_bufs[2];
    int num_active;

    sc_array_init (&recv_bufs[0], sizeof (p4est_quadrant_t));
    sc_array_init (&recv_bufs[1], sizeof (p4est_quadrant_t));
    sc_array_init (&send_bufs[0], sizeof (p4est_quadrant_t));
    sc_array_init (&send_bufs[1], sizeof (p4est_quadrant_t));

    bit[0] = bit[1] = -1;
    for (i = 0; i < 2; i++) {
      int shift = SC_LOG2_64 (procrange[i][1] - procrange[i][0]);
      if (shift >= 0) {
        bit[i] = 1 << shift;
      }
      else {
        bit[i] = 0;
      }
      state[i] = P4EST_BALSORT_ISEND;
    }

    num_active = 2;
    while (bit[0] > 0 || bit[1] > 0) {
      for (i = 0; i < 2; i++) {
        int flg;
        int rcount;
        MPI_Status status;

        if (!bit[i]) {
          num_active = 1;
          continue;
        }

        switch (state[i]) {
        case P4EST_BALSORT_ISEND:
          /* all receives have been merged, all sends have completed:
           * initiate another send and receive */
          if (num_sort == 1 && bit[i^1] > bit[i]) {
            num_active = 1;
            state[i] = P4EST_BALSORT_SYNC;
          }
          proc[i] = i ? rank + bit[i] : rank - bit[i];
          if (proc[i] >= procrange[i][0] && proc[i] <= procrange[i][1]) {
            /* prepare the send buffer */
            ssize_t sz;
            p4est_quadrant_t nextq = p4est->global_first_position[i ? proc[i] : proc[i] + 1];

            while (nextq.level > 0 && p4est_quadrant_child_id (&nextq) == 0) {
              nextq.level--;
            }
            sz = p4est_find_lower_bound (&sort_bufs[i], &nextq, 0);
            if (sz < 0) {
              sz = sort_bufs[i].elem_count;
            }
            if (!i) {
              sc_array_push_cut (&sort_bufs[i], 0, (size_t) sz, &send_bufs[i]);
            }
            else {
              sc_array_push_cut (&sort_bufs[i], (size_t) sz, sort_bufs[i].elem_count, &send_bufs[i]);
            }
            mpiret = MPI_Isend (send_bufs[i].array, send_bufs[i].elem_count * sizeof (p4est_quadrant_t), MPI_BYTE,
                                proc[i], P4EST_COMM_BALANCE_SORT_SORT, p4est->mpicomm, &sendreq[i]);
            SC_CHECK_MPI (mpiret);
            if (p4est->inspect && p4est->inspect->stats) {
              sc_statistics_accumulate (p4est->inspect->stats, "balance sort send", (double) send_bufs[i].elem_count);
            }
            state[i] = P4EST_BALSORT_IPROBE;
          }
          else {
            bit[i] >>= 1;
          }
          break;
        case P4EST_BALSORT_SYNC:
          if (bit[i ^ 1] <= bit[i]) {
            num_active = 2;
            state[i] = P4EST_BALSORT_ISEND;
          }
          break;
        case P4EST_BALSORT_IPROBE:
          /* a send has been initiated, probe for a matching receive */
          P4EST_ASSERT (proc[i] >= procrange[i][0] && proc[i] <= procrange[i][1]);
          if (num_active == 1) {
            mpiret = MPI_Probe (proc[i], P4EST_COMM_BALANCE_SORT_SORT, p4est->mpicomm,
                                &status);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret = MPI_Iprobe (proc[i], P4EST_COMM_BALANCE_SORT_SORT, p4est->mpicomm,
                                 &flg, &status);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            p4est_quadrant_t *dest;

            mpiret = MPI_Get_count (&status, MPI_BYTE, &rcount);
            SC_CHECK_MPI (mpiret);
            dest = (p4est_quadrant_t *) sc_array_push_count (&recv_bufs[i], (size_t) (rcount / sizeof (p4est_quadrant_t)));
            mpiret = MPI_Irecv (dest, rcount, MPI_BYTE, proc[i], P4EST_COMM_BALANCE_SORT_SORT,
                                p4est->mpicomm, &recvreq[i]);
            SC_CHECK_MPI (mpiret);
            if (p4est->inspect && p4est->inspect->stats) {
              sc_statistics_accumulate (p4est->inspect->stats, "balance sort recv", (double) (rcount / sizeof (p4est_quadrant_t)));
            }
            state[i] = P4EST_BALSORT_IRECV;
          }
          break;
        case P4EST_BALSORT_IRECV:
          /* a recv request has been initiated, wait for it to receive */
          if (num_active == 1) {
            mpiret = MPI_Wait (&recvreq[i], MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret = MPI_Test (&recvreq[i], &flg, MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            /* merge sort in the receive */
            if (num_sort == 2) {
              /* everything that I have received is for me */
              p4est_quadrant_array_merge_reduce (&sort_bufs[i], &recv_bufs[i]);
            }
            else {
              /* what I have received goes to the opposite side */
              p4est_quadrant_array_merge_reduce (&sort_bufs[i^1], &recv_bufs[i]);
            }
            sc_array_truncate (&recv_bufs[i]);
            state[i] = P4EST_BALSORT_WAIT;
          }
          break;
        case P4EST_BALSORT_WAIT:
          /* wait for a send to finish */
          if (num_active == 1) {
            mpiret = MPI_Wait (&sendreq[i], MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret = MPI_Test (&sendreq[i], &flg, MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            bit[i] >>= 1;
            sc_array_truncate (&send_bufs[i]);
            state[i] = P4EST_BALSORT_ISEND;
          }
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
    }

    if (num_sort == 1) {
      p4est_quadrant_array_merge_reduce (&sort_bufs[0], &sort_bufs[1]);
      sc_array_truncate (&sort_bufs[1]);
    }

    sc_array_reset (&recv_bufs[1]);
    sc_array_reset (&recv_bufs[0]);
    sc_array_reset (&send_bufs[1]);
    sc_array_reset (&send_bufs[0]);
  }
#endif
  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort_merge_seeds (p4est_t *p4est, p4est_connect_type_t btype,
                                int flt, int llt, int num_trees,
                                int num_sort, int *minlevel,
                                p4est_quadrant_t *desc,
                                sc_array_t *tree_bufs, sc_array_t *sort_bufs)
{
  sc_flopinfo_t snap;

  P4EST_FUNC_SNAP (p4est, &snap);
#ifdef P4EST_ENABLE_MPI
  if (num_sort > 0) {
    int i;

    /* combine sort results back into the original trees */
    for (i = 0; i < num_sort; i++) {
      p4est_quadrant_t anc;
      size_t nq, s;
      sc_array_t seeds;
      sc_array_t newquads;

      p4est_quadrant_ancestor (&desc[i], minlevel[i], &anc);

      sc_array_init (&seeds, sizeof (p4est_quadrant_t));
      sc_array_init (&newquads, sizeof (p4est_quadrant_t));
      nq = sort_bufs[i].elem_count;
      for (s = 0; s < nq; s++) {
        p4est_quadrant_t *q = (p4est_quadrant_t *) sc_array_index (&sort_bufs[i], s);
        p4est_quadrant_t *prev = s ? (p4est_quadrant_t *) sc_array_index (&sort_bufs[i], s - 1) : NULL;
        p4est_quadrant_t *next = (s < nq - 1) ? (p4est_quadrant_t *) sc_array_index (&sort_bufs[i], s + 1) : NULL;
        p4est_quadrant_t lb, ub;
        int32_t touch;
        int32_t mask = 1;

        P4EST_ASSERT (p4est_quadrant_is_ancestor (&anc, q));
        if (prev) {
          p4est_quadrant_t a, qstop;
          int qid, pid, idxor, dir, qstopid;

          p4est_nearest_common_ancestor (q, prev, &a);
          P4EST_ASSERT (a.level < q->level && a.level < prev->level);
          qid = p4est_quadrant_ancestor_id (q, a.level + 1);
          pid = p4est_quadrant_ancestor_id (prev, a.level + 1);
          P4EST_ASSERT (qid > pid);
          idxor = qid ^ pid;
          dir = SC_LOG2_8 (idxor);
          P4EST_ASSERT (dir >= 0);
          qstopid = (qid >> dir) << dir;
          p4est_quadrant_child (&a, &qstop, qstopid);
          p4est_quadrant_first_descendant (&qstop, &lb, P4EST_QMAXLEVEL);
        }
        if (next) {
          p4est_quadrant_t a, qstop;
          int qid, nid, idxor, dir, qstopid;

          p4est_nearest_common_ancestor (q, next, &a);
          P4EST_ASSERT (a.level < q->level && a.level < next->level);
          qid = p4est_quadrant_ancestor_id (q, a.level + 1);
          nid = p4est_quadrant_ancestor_id (next, a.level + 1);
          P4EST_ASSERT (nid > qid);
          idxor = qid ^ nid;
          dir = SC_LOG2_8 (idxor);
          P4EST_ASSERT (dir >= 0);
          qstopid = ((nid >> dir) << dir) - 1;
          p4est_quadrant_child (&a, &qstop, qstopid);
          p4est_quadrant_last_descendant (&qstop, &ub, P4EST_QMAXLEVEL);
        }
        touch = p4est_find_range_boundaries (prev ? &lb : NULL, next ? &ub : NULL,
                                             anc.level, NULL,
#ifdef P4_TO_P8
                                             NULL,
#endif
                                             NULL);
        {
          int f;

          for (f = 0; f < P4EST_FACES; f++, mask <<= 1) {
            p4est_quadrant_t n;
            int split;

            if (!(touch & mask)) {
              continue;
            }
            p4est_quadrant_face_neighbor (&anc, f, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
#ifdef P4_TO_P8
        {
          int e;

          for (e = 0; e < P8EST_EDGES; e++, mask <<= 1) {
            p4est_quadrant_t n;
            int split;

            if (!(touch & mask)) {
              continue;
            }
            p8est_quadrant_edge_neighbor (&anc, e, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
#endif
        {
          int c;

          for (c = 0; c < P4EST_CHILDREN; c++, mask <<= 1) {
            p4est_quadrant_t n;
            int split;

            if (!(touch & mask)) {
              continue;
            }
            p4est_quadrant_corner_neighbor (&anc, c, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
      }
      sc_array_sort (&newquads, p4est_quadrant_compare);
      p4est_quadrant_array_reduce (&newquads, NULL, NULL);
      p4est_quadrant_array_merge_reduce (&sort_bufs[i], &newquads);
      p4est_quadrant_array_merge_reduce (&tree_bufs[desc[i].p.which_tree - flt], &sort_bufs[i]);
      sc_array_reset (&seeds);
      sc_array_reset (&newquads);
    }

  }
#endif
  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort_neigh (p4est_t *p4est, int flt, int llt, int num_trees,
                          int *proc_hash, sc_array_t *procs,
                          p4est_tree_neigh_info_t *tinfo,
                          sc_array_t *tree_bufs,
                          sc_array_t *neigh_bufs)
{
  p4est_topidx_t t;
  sc_flopinfo_t snap;
  size_t nneigh = procs->elem_count;
  int rank = p4est->mpirank;

  P4EST_FUNC_SNAP (p4est, &snap);
  /* ghost neighbor balance */
  for (t = 0; t < num_trees; t++) {
    int i, j, k, l;
#ifdef P4_TO_P8
    int kstart = 0, kend = 3;
#else
    int kstart = 1, kend = 2;
#endif
    sc_array_t *tquads = &tree_bufs[t];
    sc_array_t orig_view;
    sc_array_t tform_quads;
    p4est_qcoord_t H = P4EST_ROOT_LEN;
    p4est_tree_neigh_info_t *info = &tinfo[t];
    p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, t + flt);

    sc_array_init (&tform_quads, sizeof (p4est_quadrant_t));
    for (k = kstart, l = 0; k < kend; k++) {
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++, l++) {
          p4est_quadrant_t insul, ub;
          size_t sz;
          size_t qstart, qend;

          insul.x = (i - 1) * H;
          insul.y = (j - 1) * H;
#ifdef P4_TO_P8
          insul.z = (k - 1) * H;
#endif
          insul.level = 0;
          p4est_quadrant_last_descendant (&insul, &ub, P4EST_QMAXLEVEL);
          sz = p4est_find_lower_bound (tquads, &insul, 0);
          if (sz < 0) {
            sz = tquads->elem_count;
          }
          qstart = (size_t) sz;
          sz = p4est_find_higher_bound (tquads, &ub, 0);
          if (sz < 0) {
            sz = 0;
          }
          else {
            sz++;
          }
          qend = (size_t) sz;
          if (qend <= qstart) {
            continue;
          }
          sc_array_init_view (&orig_view, tquads, qstart, qend - qstart);
          if (l == P4EST_INSUL / 2) {
            p4est_balance_sort_divide (p4est, procs, neigh_bufs, t + flt, &orig_view, 1);
          }
          else {
            size_t z, sq;
            size_t no = orig_view.elem_count;

            for (z = info->offset[l]; z < info->offset[l + 1]; z++) {
              p4est_tree_neigh_t *tn = (p4est_tree_neigh_t *) sc_array_index (&info->tnarray, z);
              p4est_quadrant_t *oq = (p4est_quadrant_t *) sc_array_index (&orig_view, 0);
              p4est_quadrant_t *trq;

              trq = (p4est_quadrant_t *) sc_array_push_count (&tform_quads, no);

              for (sq = 0; sq < no; sq++) {
                p4est_quadrant_t temp;
                p4est_quadrant_utransform (&oq[sq], &temp, &(tn->u[0]), 0);

                p4est_quadrant_sibling (&temp, &trq[sq], 0);
                P4EST_ASSERT (p4est_quadrant_is_valid (&trq[sq]));
              }
              sc_array_sort (&tform_quads, p4est_quadrant_compare);
              p4est_balance_sort_divide (p4est, procs, neigh_bufs, tn->nt, &tform_quads, 0);
              sc_array_truncate (&tform_quads);
            }
          }
        }
      }
    }
    {
      p4est_quadrant_t fd = tree->first_desc;
      ssize_t sz;
      size_t qstart, qend;

      while (fd.level > 0 && !p4est_quadrant_child_id (&fd)) {
        fd.level--;
      }
      sz = p4est_find_lower_bound (tquads, &fd, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qstart = sz;
      sz = p4est_find_higher_bound (tquads, &tree->last_desc, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      qend = sz;
      if (qend > qstart) {
        sc_array_crop (tquads, qstart, qend);
      }
      else {
        sc_array_truncate (tquads);
      }
    }
    sc_array_reset (&tform_quads);
  }
  /* exchange, merge and complete */
  if (num_trees) {
#ifdef P4EST_ENABLE_MPI
    size_t n;
    MPI_Request *req;
    sc_array_t *recv_buf = NULL;
    int me = proc_hash[rank];
    int mpiret;

    P4EST_ASSERT (me > 0 && me <= nneigh);

    req = P4EST_ALLOC (MPI_Request, nneigh);

    req[me - 1] = MPI_REQUEST_NULL;

    recv_buf = &neigh_bufs[me - 1];
    for (n = 0; n < nneigh; n++) {
      int p = *((int *) sc_array_index (procs, n));
      sc_array_t *buf = &neigh_bufs[n];

      if (p == rank) {
        P4EST_ASSERT (n == me - 1);
        continue;
      }
      mpiret = MPI_Isend (buf->array, buf->elem_count * sizeof (p4est_quadrant_t), MPI_BYTE,
                          p, P4EST_COMM_BALANCE_SORT_NEIGH, p4est->mpicomm, &req[n]);
      SC_CHECK_MPI (mpiret);
    }
    for (n = 0; n < nneigh - 1; n++) {
      MPI_Status status;
      int rcount;
      int p;
      p4est_quadrant_t *q;

      mpiret = MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_BALANCE_SORT_NEIGH,
                          p4est->mpicomm, &status);
      SC_CHECK_MPI (mpiret);
      p = status.MPI_SOURCE;
      P4EST_ASSERT (proc_hash[p] > 0);
      mpiret = MPI_Get_count (&status, MPI_BYTE, &rcount);
      P4EST_ASSERT ((rcount % sizeof (p4est_quadrant_t)) == 0);
      q = (p4est_quadrant_t *) sc_array_push_count (recv_buf, rcount / sizeof (p4est_quadrant_t));
      mpiret = MPI_Recv (q, rcount, MPI_BYTE, p, P4EST_COMM_BALANCE_SORT_NEIGH,
                         p4est->mpicomm, MPI_STATUS_IGNORE);
      SC_CHECK_MPI (mpiret);
    }
    mpiret = MPI_Waitall ((int) nneigh, req, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_FREE (req);
#endif

    sc_array_sort (recv_buf, p4est_quadrant_compare_piggy);
#ifdef P4EST_ENABLE_DEBUG
    if (recv_buf->elem_count) {
      p4est_quadrant_t *first = p4est_quadrant_array_index (recv_buf, 0);
      p4est_quadrant_t *last = p4est_quadrant_array_index (recv_buf, recv_buf->elem_count - 1);
      p4est_quadrant_t *gf = &p4est->global_first_position[rank];
      p4est_quadrant_t *gn = &p4est->global_first_position[rank + 1];

      P4EST_ASSERT (first->p.which_tree >= gf->p.which_tree);
      P4EST_ASSERT (last->p.which_tree <= gn->p.which_tree);
      if (first->p.which_tree == gf->p.which_tree) {
        P4EST_ASSERT (p4est_quadrant_compare (gf, first) < 0
                      ||
                      (gf->x == first->x && gf->y == first->y &&
#ifdef P4_TO_P8
                       gf->z == first->z &&
#endif
                       1));

      }
      if (last->p.which_tree == gn->p.which_tree) {
        P4EST_ASSERT (p4est_quadrant_compare (last, gn) < 0);
      }
    }
#endif
    {
      size_t s, n;

      n = recv_buf->elem_count;
      for (s = 0; s < n;) {
        size_t z;
        p4est_quadrant_t *q = p4est_quadrant_array_index (recv_buf, s);
        p4est_topidx_t t = q->p.which_tree;
        sc_array_t view;

        for (z = s + 1; z < n; z++) {
          p4est_quadrant_t *r = p4est_quadrant_array_index (recv_buf, z);

          if (r->p.which_tree != t) {
            break;
          }
        }
        sc_array_init_view (&view, recv_buf, s, z - s);
        p4est_quadrant_array_reduce (&view, NULL, NULL);
        p4est_quadrant_array_merge_reduce (&tree_bufs[t - flt], &view);
        s = z;
      }
    }
  }
  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort_complete (p4est_t *p4est, int flt, int llt, int num_trees,
                             sc_array_t *tree_bufs,
                             p4est_init_t init_fn,
                             p4est_replace_t replace_fn)
{
  p4est_topidx_t t;
  sc_flopinfo_t snap;

  P4EST_FUNC_SNAP (p4est, &snap);

  if (num_trees) {
    sc_array_t newcomplete;

    sc_array_init (&newcomplete, sizeof (p4est_quadrant_t));
    p4est->local_num_quadrants = 0;
    for (t = flt; t <= llt; t++) {
      p4est_tree_t *tree = p4est_tree_array_index (p4est->trees, t);
      p4est_quadrant_t dom;

      p4est_nearest_common_ancestor (&tree->first_desc, &tree->last_desc, &dom);
      (void) sc_array_push_count (&newcomplete, tree_bufs[t-flt].elem_count);
      sc_array_truncate (&newcomplete);
      p4est_complete_kernel (&tree_bufs[t - flt], &dom, &tree->first_desc, &tree->last_desc, &newcomplete);
      p4est_subtree_replace (p4est, t, &newcomplete, init_fn, replace_fn);
      p4est->local_num_quadrants += newcomplete.elem_count;
      sc_array_truncate (&newcomplete);
    }
    sc_array_reset (&newcomplete);
  }

  P4EST_FUNC_SHOT (p4est, &snap);
}

static void
p4est_balance_sort (p4est_t *p4est, p4est_connect_type_t btype,
                    p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_topidx_t t, num_trees;
  p4est_topidx_t flt = p4est->first_local_tree;
  p4est_topidx_t llt = p4est->last_local_tree;
  p4est_tree_neigh_info_t *tinfo = NULL;
  int num_sort = -1;
  sc_flopinfo_t snap;
  sc_array_t         *procs;
  int *proc_hash;
  int minlevel[2];
  int procrange[2][2];
  p4est_quadrant_t desc[2];
  size_t s, nneigh;
  int rank = p4est->mpirank;
  sc_array_t *neigh_bufs;
  sc_array_t *tree_bufs;
  sc_array_t sort_bufs[2];

  P4EST_FUNC_SNAP (p4est, &snap);

  /* first figure out the parallel neighborhood */
  num_trees = SC_MAX (0, llt + 1 - flt);
  if (!num_trees) {
    p4est_quadrant_t q;

    q = p4est->global_first_position[rank];
    P4EST_ASSERT (p4est_quadrant_is_equal_piggy (&q, &p4est->global_first_position[rank + 1]));
    if (q.x != 0 || q.y != 0 ||
#ifdef P4_TO_P8
        q.z != 0 ||
#endif
        0) {
      flt = llt = q.p.which_tree;
    }
  }
  proc_hash = P4EST_ALLOC_ZERO (int, p4est->mpisize);
  procs = sc_array_new_size (sizeof (int), P4EST_INSUL);
  tinfo = P4EST_ALLOC (p4est_tree_neigh_info_t, llt + 1 - flt);
  p4est_balance_sort_compute_pattern (p4est, flt, llt, num_trees,
                                      tinfo, proc_hash, procs,
                                      procrange, minlevel, &num_sort, desc);
  nneigh = procs->elem_count;
  neigh_bufs = P4EST_ALLOC (sc_array_t, nneigh);
  for (s = 0; s < nneigh; s++) {
    sc_array_init (&neigh_bufs[s], sizeof (p4est_quadrant_t));
  }
  tree_bufs = P4EST_ALLOC (sc_array_t, num_trees);
  for (t = 0; t < num_trees; t++) {
    sc_array_init (&tree_bufs[t], sizeof (p4est_quadrant_t));
  }

  sc_array_init (&sort_bufs[0], sizeof (p4est_quadrant_t));
  sc_array_init (&sort_bufs[1], sizeof (p4est_quadrant_t));

  p4est_balance_sort_local (p4est, btype, flt, llt, num_trees,
                            num_sort, minlevel, desc, tree_bufs,
                            sort_bufs);

  p4est_balance_sort_sort (p4est, flt, llt, num_trees, num_sort,
                           procrange, sort_bufs);

  p4est_balance_sort_merge_seeds (p4est, btype, flt, llt, num_trees,
                                  num_sort, minlevel, desc,
                                  tree_bufs, sort_bufs);

  sc_array_reset (&sort_bufs[1]);
  sc_array_reset (&sort_bufs[0]);
  p4est_balance_sort_neigh (p4est, flt, llt, num_trees, proc_hash, procs, tinfo,
                            tree_bufs, neigh_bufs);

  p4est_balance_sort_complete (p4est, flt, llt, num_trees, tree_bufs,
                               init_fn, replace_fn);

  for (t = 0; t < num_trees; t++) {
    sc_array_reset (&tree_bufs[t]);
  }
  P4EST_FREE (tree_bufs);
  for (s = 0; s < nneigh; s++) {
    sc_array_reset (&neigh_bufs[s]);
  }
  P4EST_FREE (neigh_bufs);
  sc_array_destroy (procs);
  for (t = flt; t <= llt; t++) {
    sc_array_reset (&(tinfo[t - flt].tnarray));
  }
  P4EST_FREE (tinfo);
  P4EST_FREE (proc_hash);
  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_balance_ext (p4est_t * p4est, p4est_connect_type_t btype,
                   p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_gloidx_t      old_gnq;
  const int8_t       *pre_adapt_flags = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_balance %s with %lld total quadrants\n",
                            p4est_connect_type_string (btype),
                            (long long) p4est->global_num_quadrants);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not decrease */
  old_gnq = p4est->global_num_quadrants;
  if (p4est->inspect && p4est->inspect->balance_sort) {
    p4est_balance_sort (p4est, btype, NULL, NULL);
  }
  else {
    p4est_balance_ext_dirty (p4est, btype, init_fn, replace_fn);
  }
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  if (p4est->inspect) {
    pre_adapt_flags = p4est->inspect->pre_adapt_flags;
  }
  P4EST_ASSERT (pre_adapt_flags || p4est->global_num_quadrants >= old_gnq);
  if (pre_adapt_flags || old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est_is_balanced (p4est, btype));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_balance with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_balance_ext_dirty (p4est_t * p4est, p4est_connect_type_t btype,
                         p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  const int           rank = p4est->mpirank;
  const int           num_procs = p4est->mpisize;
  int                 j, k, l, m, which;
#ifdef P4EST_ENABLE_DEBUG
  int                 facewhich;
#endif
  int                 face;
  int                 first_peer, last_peer;
  int                 quad_contact[P4EST_FACES];
  int                 any_boundary, tree_contact[P4EST_INSUL];
  int                 tree_fully_owned, full_tree[2];
  int8_t             *tree_flags;
  size_t              zz, treecount, ctree;
  size_t              localcount;
  size_t              qcount, qbytes;
  size_t              all_incount, all_outcount;
  p4est_qcoord_t      qh;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_topidx_t      qtree, nt;
  p4est_topidx_t      first_tree, last_tree;
  p4est_locidx_t      skipped;
  p4est_balance_peer_t *peers, *peer;
  p4est_tree_t       *tree;
  p4est_quadrant_t    mylow, nextlow;
  p4est_quadrant_t    tosend, insulq, tempq;
  p4est_quadrant_t   *q, *s;
  p4est_connectivity_t *conn = p4est->connectivity;
  sc_array_t         *qarray, *tquadrants;
  sc_array_t         *borders;
  const int8_t       *pre_adapt_flags = NULL;
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  int                 ftransform[P4EST_FTRANSFORM];
  int                 face_axis[3];     /* 3 not P4EST_DIM */
  int                 contact_face_only, contact_edge_only;
#ifdef P4_TO_P8
  int                 edge;
  size_t              etree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
#endif
  int                 corner;
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;
#ifdef P4EST_ENABLE_MPI
#ifdef P4EST_ENABLE_DEBUG
  unsigned            checksum;
  sc_array_t          checkarray;
  p4est_gloidx_t      ltotal[2], gtotal[2];
#endif /* P4EST_ENABLE_DEBUG */
  int                 i;
  int                 mpiret, rcount;
  int                 request_first_send, request_first_recv, outcount;
  int                 request_second_send, total_send_count, total_recv_count;
  int                 send_zero[2], send_load[2];
  int                 recv_zero[2], recv_load[2];
  int                *wait_indices;
  int                *receiver_ranks, *sender_ranks;
  int                 num_receivers, num_senders;
  MPI_Request        *send_requests_first_load;
  MPI_Request        *send_requests_second_load;
  MPI_Request        *recv_requests_first_load;
  MPI_Request        *recv_requests_second_load;
  MPI_Status         *recv_statuses, *jstatus;
  sc_array_t         *receivers, *senders;
  sc_array_t         *in_counts, *out_counts;
  int                *icounts;
  sc_flopinfo_t       snap;
#endif /* P4EST_ENABLE_MPI */
  sc_notify_t        *notify = NULL;
  int                 own_notify = 0;

  P4EST_FUNC_SNAP (p4est, &snap);

#ifndef P4_TO_P8
  P4EST_ASSERT (btype == P4EST_CONNECT_FACE || btype == P4EST_CONNECT_CORNER);
#else
  P4EST_ASSERT (btype == P8EST_CONNECT_FACE || btype == P8EST_CONNECT_EDGE ||
                btype == P8EST_CONNECT_CORNER);
#endif

#ifdef P4EST_ENABLE_DEBUG
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
#endif

  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&tosend);
  P4EST_QUADRANT_INIT (&insulq);
  P4EST_QUADRANT_INIT (&tempq);

  /* tree status flags (max 8 per tree) */
  tree_flags = P4EST_ALLOC (int8_t, conn->num_trees);
  for (nt = 0; nt < conn->num_trees; ++nt) {
    tree_flags[nt] = 0x00;
  }

  localcount = (size_t) (p4est->last_local_tree + 1 -
                         p4est->first_local_tree);
  borders = sc_array_new_size (sizeof (sc_array_t), localcount);
  for (zz = 0; zz < localcount; zz++) {
    qarray = (sc_array_t *) sc_array_index (borders, zz);
    sc_array_init (qarray, sizeof (p4est_quadrant_t));
  }

#ifdef P4EST_ENABLE_MPI
  send_requests_first_load = P4EST_ALLOC (MPI_Request, 4 * num_procs);
  recv_requests_first_load = send_requests_first_load + 1 * num_procs;
  send_requests_second_load = send_requests_first_load + 2 * num_procs;
  recv_requests_second_load = send_requests_first_load + 3 * num_procs;
  recv_statuses = P4EST_ALLOC (MPI_Status, num_procs);
  for (j = 0; j < num_procs; ++j) {
    send_requests_first_load[j] = MPI_REQUEST_NULL;
    recv_requests_first_load[j] = MPI_REQUEST_NULL;
    send_requests_second_load[j] = MPI_REQUEST_NULL;
    recv_requests_second_load[j] = MPI_REQUEST_NULL;
  }
  wait_indices = P4EST_ALLOC (int, num_procs);
#ifdef P4EST_ENABLE_DEBUG
  sc_array_init (&checkarray, 4);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* allocate per peer storage and initialize requests */
  peers = P4EST_ALLOC (p4est_balance_peer_t, num_procs);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    sc_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_second, sizeof (p4est_quadrant_t));
    peer->send_first_count = peer->send_second_count = 0;
    peer->recv_first_count = peer->recv_second_count = 0;
    peer->have_first_count = peer->have_first_load = 0;
    peer->have_second_count = peer->have_second_load = 0;
  }
#ifdef P4_TO_P8
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  /* compute first quadrant on finest level */
  mylow.x = p4est->global_first_position[rank].x;
  mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
  mylow.z = p4est->global_first_position[rank].z;
#endif
  mylow.level = P4EST_QMAXLEVEL;

  /* and the first finest quadrant of the next processor */
  nextlow.x = p4est->global_first_position[rank + 1].x;
  nextlow.y = p4est->global_first_position[rank + 1].y;
#ifdef P4_TO_P8
  nextlow.z = p4est->global_first_position[rank + 1].z;
#endif
  nextlow.level = P4EST_QMAXLEVEL;

  /* start balance_A timing */
  if (p4est->inspect != NULL) {
    p4est->inspect->balance_A = -sc_MPI_Wtime ();
    p4est->inspect->balance_A_count_in = 0;
    p4est->inspect->balance_A_count_out = 0;
    p4est->inspect->use_B = 0;
    pre_adapt_flags = p4est->inspect->pre_adapt_flags;
  }

  /* loop over all local trees to assemble first send list */
  first_tree = p4est->first_local_tree;
  last_tree = p4est->last_local_tree;
  first_peer = rank;
  last_peer = rank;
  all_incount = 0;
  skipped = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    const int8_t       *this_pre_adapt_flags = NULL;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);
    tree_fully_owned = full_tree[0] && full_tree[1];
    any_boundary = 0;
    for (i = 0; i < P4EST_INSUL; ++i) {
      any_boundary |= tree_contact[i];
    }
    if (any_boundary) {
      tree_flags[nt] |= any_face_flag;
    }
    tree = p4est_tree_array_index (p4est->trees, nt);
    tquadrants = &tree->quadrants;
    if (pre_adapt_flags) {
      this_pre_adapt_flags = &pre_adapt_flags[all_incount];
    }
    all_incount += tquadrants->elem_count;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into balance tree %lld with %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);

    /* local balance first pass */
    p4est_balance_subtree_ext (p4est, btype, nt, this_pre_adapt_flags,
                               init_fn, replace_fn);
    treecount = tquadrants->elem_count;
    P4EST_VERBOSEF ("Balance tree %lld A %llu\n",
                    (long long) nt, (unsigned long long) treecount);

    /* check if this tree is not shared with other processors */
    if (tree_fully_owned) {
      /* all quadrants in this tree are owned by me */
      tree_flags[nt] |= fully_owned_flag;
      if (!any_boundary) {
        /* this tree is isolated, no balance between trees */
        continue;
      }
    }

    if (borders != NULL) {
      qarray = (sc_array_t *) sc_array_index (borders,
                                              (size_t) (nt - first_tree));
    }
    else {
      qarray = NULL;
    }

    /* identify boundary quadrants and prepare them to be sent */
    for (zz = 0; zz < treecount; ++zz) {
      /* this quadrant may be on the boundary with a range of processors */
      q = p4est_quadrant_array_index (tquadrants, zz);
      qh = P4EST_QUADRANT_LEN (q->level);
      if (p4est_comm_neighborhood_owned (p4est, nt,
                                         full_tree, tree_contact, q)) {
        /* this quadrant's 3x3 neighborhood is owned by this processor */
        ++skipped;
        continue;
      }

      if (qarray != NULL) {
        s = (p4est_quadrant_t *) sc_array_push (qarray);
        *s = *q;
      }

#ifdef P4_TO_P8
      for (m = 0; m < 3; ++m) {
#if 0
      }
#endif
#else
      m = 0;
#endif
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) {
          which = m * 9 + k * 3 + l;    /* 2D: 0..8, 3D: 0..26 */
          /* exclude myself from the queries */
          if (which == P4EST_INSUL / 2) {
            continue;
          }
          /* may modify insulq below, never modify q itself! */
          insulq = *q;
          insulq.x += (l - 1) * qh;
          insulq.y += (k - 1) * qh;
#ifdef P4_TO_P8
          insulq.z += (m - 1) * qh;
#endif
          /* check boundary status of insulation quadrant */
          quad_contact[0] = (insulq.x < 0);
          quad_contact[1] = (insulq.x >= rh);
          face_axis[0] = quad_contact[0] || quad_contact[1];
          quad_contact[2] = (insulq.y < 0);
          quad_contact[3] = (insulq.y >= rh);
          face_axis[1] = quad_contact[2] || quad_contact[3];
#ifndef P4_TO_P8
          face_axis[2] = 0;
#else
          quad_contact[4] = (insulq.z < 0);
          quad_contact[5] = (insulq.z >= rh);
          face_axis[2] = quad_contact[4] || quad_contact[5];
          edge = -1;
#endif
          contact_edge_only = contact_face_only = 0;
          face = -1;
          if (face_axis[0] || face_axis[1] || face_axis[2]) {
            /* this quadrant is relevant for inter-tree balancing */
            if (!face_axis[1] && !face_axis[2]) {
              contact_face_only = 1;
              face = 0 + quad_contact[1];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + (2 * quad_contact[1] - 1);
#endif
            }
            else if (!face_axis[0] && !face_axis[2]) {
              contact_face_only = 1;
              face = 2 + quad_contact[3];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + 3 * (2 * quad_contact[3] - 1);
#endif
            }
#ifdef P4_TO_P8
            else if (!face_axis[0] && !face_axis[1]) {
              contact_face_only = 1;
              face = 4 + quad_contact[5];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + 9 * (2 * quad_contact[5] - 1);
#endif
            }
            else if (!face_axis[0]) {
              contact_edge_only = 1;
              edge = 0 + 2 * quad_contact[5] + quad_contact[3];
            }
            else if (!face_axis[1]) {
              contact_edge_only = 1;
              edge = 4 + 2 * quad_contact[5] + quad_contact[1];
            }
            else if (!face_axis[2]) {
              contact_edge_only = 1;
              edge = 8 + 2 * quad_contact[3] + quad_contact[1];
            }
#endif
            if (contact_face_only) {
              /* square contact across a face */
              P4EST_ASSERT (!contact_edge_only);
              P4EST_ASSERT (face >= 0 && face < P4EST_FACES);
              P4EST_ASSERT (quad_contact[face]);
              qtree = p4est_find_face_transform (conn, nt, face, ftransform);
              if (qtree >= 0) {
                P4EST_ASSERT (tree_contact[facewhich]);
                p4est_quadrant_transform_face (q, &tosend, ftransform);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = face;
                p4est_quadrant_transform_face (&insulq, &tempq, ftransform);
                p4est_balance_schedule (p4est, peers, qtree, 1,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
              else {
                /* goes across a face with no neighbor */
                P4EST_ASSERT (!tree_contact[facewhich]);
              }
            }
#ifdef P4_TO_P8
            else if (contact_edge_only) {
              /* this quadrant crosses an edge */
              P4EST_ASSERT (!contact_face_only);
              P4EST_ASSERT (edge >= 0 && edge < P8EST_EDGES);
              p8est_find_edge_transform (conn, nt, edge, &ei);
              for (etree = 0; etree < eta->elem_count; ++etree) {
                et = p8est_edge_array_index (eta, etree);
                p8est_quadrant_transform_edge (q, &tosend, &ei, et, 0);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = edge;
                p8est_quadrant_transform_edge (&insulq, &tempq, &ei, et, 1);
                p4est_balance_schedule (p4est, peers, et->ntree, 1,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
            }
#endif
            else {
              /* this quadrant crosses a corner */
              P4EST_ASSERT (face_axis[0] && face_axis[1]);
              corner = quad_contact[1] + 2 * quad_contact[3];
#ifdef P4_TO_P8
              P4EST_ASSERT (face_axis[2]);
              corner += 4 * quad_contact[5];
#endif
              P4EST_ASSERT (p4est_quadrant_touches_corner (q, corner, 1));
              P4EST_ASSERT (p4est_quadrant_touches_corner
                            (&insulq, corner, 0));
              p4est_find_corner_transform (conn, nt, corner, &ci);
              for (ctree = 0; ctree < cta->elem_count; ++ctree) {
                ct = p4est_corner_array_index (cta, ctree);
                tosend = *q;
                p4est_quadrant_transform_corner (&tosend, (int) ct->ncorner,
                                                 0);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = corner;
                tempq = insulq;
                p4est_quadrant_transform_corner (&tempq, (int) ct->ncorner,
                                                 1);
                p4est_balance_schedule (p4est, peers, ct->ntree, 1,
                                        &tosend, &tempq, &first_peer,
                                        &last_peer);
              }
            }
          }
          else {
            /* no inter-tree contact */
            tosend = *q;
            tosend.p.piggy2.from_tree = nt;
            tosend.pad16 = -1;
            p4est_balance_schedule (p4est, peers, nt, 0,
                                    &tosend, &insulq, &first_peer,
                                    &last_peer);
          }
        }
      }
#ifdef P4_TO_P8
#if 0
      {
#endif
      }
#endif
    }
    tquadrants = NULL;          /* safeguard */
  }

  /* end balance_A, start balance_comm */
  if (p4est->inspect != NULL) {
    p4est->inspect->balance_A += sc_MPI_Wtime ();
    p4est->inspect->balance_comm = -sc_MPI_Wtime ();
    p4est->inspect->balance_comm_sent = 0;
    p4est->inspect->balance_comm_nzpeers = 0;
    for (k = 0; k < 2; ++k) {
      p4est->inspect->balance_load_sends[k] = 0;
      p4est->inspect->balance_load_receives[k] = 0;
      p4est->inspect->balance_zero_sends[k] = 0;
      p4est->inspect->balance_zero_receives[k] = 0;
    }
    notify = p4est->inspect->notify;
  }
  if (!notify) {
    notify = sc_notify_new (p4est->mpicomm);
    own_notify = 1;
  }

#ifdef P4EST_ENABLE_MPI
  /* encode and distribute the asymmetric communication pattern */
  receivers = sc_array_new (sizeof (int));
  senders = sc_array_new (sizeof (int));
  in_counts = sc_array_new (sizeof (int));
  out_counts = sc_array_new (sizeof (int));
#ifdef P4EST_ENABLE_DEBUG
  for (j = 0; j < first_peer; j++) {
    P4EST_ASSERT (peers[j].send_first.elem_count == 0);
  }
#endif
  for (num_receivers = 0, j = first_peer; j <= last_peer; j++) {
    if (j == rank) {
      continue;
    }
    if (peers[j].send_first.elem_count) {
      num_receivers++;
    }
  }
#ifdef P4EST_ENABLE_DEBUG
  for (j = last_peer + 1; j < num_procs; j++) {
    P4EST_ASSERT (peers[j].send_first.elem_count == 0);
  }
#endif
  sc_array_resize (receivers, (size_t) num_receivers);
  sc_array_resize (in_counts, (size_t) num_receivers);
  receiver_ranks = (int *) receivers->array;
  icounts = (int *) in_counts->array;
  for (num_receivers = 0, j = first_peer; j <= last_peer; j++) {
    if (j == rank) {
      continue;
    }
    if (peers[j].send_first.elem_count) {
      receiver_ranks[num_receivers] = j;
      icounts[num_receivers++] = (int) peers[j].send_first.elem_count;
    }
  }

  if (p4est->inspect && p4est->inspect->stats) {
    sc_notify_set_stats (notify, p4est->inspect->stats);
  }
  sc_notify_payload (receivers, senders, in_counts, out_counts, 0, notify);
  sc_array_destroy (in_counts);
  sender_ranks = (int *) senders->array;
  num_senders = (int) senders->elem_count;

  /*
   * loop over all peers and send first round of quadrants
   * for intra-tree balancing, each load is contained in one tree
   */
  total_send_count = total_recv_count = 0;
  request_first_send = 0;
  request_first_recv = 0;
  request_second_send = 0;
  send_zero[0] = send_load[0] = recv_zero[0] = recv_load[0] = 0;
  send_zero[1] = send_load[1] = recv_zero[1] = recv_load[1] = 0;

  /* Use sender_ranks array to recv from them */
  icounts = (int *) out_counts->array;
  for (k = 0; k < num_senders; k++) {
    j = sender_ranks[k];
    qcount = icounts[k];
    P4EST_ASSERT (j >= 0 && j < num_procs && j != rank);
    peer = peers + j;
    sc_array_resize (&(peer->recv_first), (size_t) qcount);
    /* first send number of quadrants to be expected */
    if (qcount > 0) {
      P4EST_LDEBUGF ("Balance A recv %llu quadrants from %d\n",
                     (unsigned long long) qcount, j);
      total_recv_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      peer->recv_first_count = qcount;
      mpiret = MPI_Irecv (peer->recv_first.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &recv_requests_first_load[j]);
      SC_CHECK_MPI (mpiret);
      ++recv_load[0];
      request_first_recv++;
    }
    else {
      ++recv_zero[0];
    }
  }
  sc_array_destroy (out_counts);
  /* Use receiver_ranks array to send to them */
  for (k = 0; k < num_receivers; ++k) {
    j = receiver_ranks[k];
    P4EST_ASSERT (j >= first_peer && j <= last_peer && j != rank);
    peer = peers + j;
    qcount = peer->send_first.elem_count;

    /* first send number of quadrants to be expected */
    if (qcount > 0) {
      P4EST_LDEBUGF ("Balance A send %llu quadrants to %d\n",
                     (unsigned long long) qcount, j);
      ++send_load[0];
    }
    else {
      ++send_zero[0];
    }
    peer->send_first_count = (int) qcount;
    /* sort and send the actual quadrants and post receive for reply */
    if (qcount > 0) {
      sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);

#ifdef P4EST_ENABLE_DEBUG
      checksum = p4est_quadrant_checksum (&peer->send_first, &checkarray, 0);
      P4EST_LDEBUGF ("Balance A send checksum 0x%08x to %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */

      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_first.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &send_requests_first_load[j]);
      SC_CHECK_MPI (mpiret);
      ++request_first_send;
    }
  }
  peer = NULL;

  /* wait for quadrant counts and post receive and send for quadrants */
  while (request_first_recv > 0) {
    mpiret = MPI_Waitsome (num_procs, recv_requests_first_load,
                           &outcount, wait_indices, recv_statuses);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      jstatus = &recv_statuses[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (j != rank && 0 <= j && j < num_procs);
      P4EST_ASSERT (jstatus->MPI_SOURCE == j);

      /* check if we are in receiving count or load */
      peer = peers + j;
      /* verify received size */
      P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_FIRST_LOAD);
      P4EST_ASSERT (peer->recv_first_count > 0);
      mpiret = MPI_Get_count (jstatus, MPI_BYTE, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount ==
                       peer->recv_first_count *
                       (int) sizeof (p4est_quadrant_t),
                       "Receive load mismatch A %d %dx%llu", rcount,
                       peer->recv_first_count,
                       (unsigned long long) sizeof (p4est_quadrant_t));

      /* received load, close this request */
      peer->have_first_load = 1;
      --request_first_recv;

#ifdef P4EST_ENABLE_DEBUG
      checksum = p4est_quadrant_checksum (&peer->recv_first, &checkarray, 0);
      P4EST_LDEBUGF ("Balance A recv checksum 0x%08x from %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */

      /* process incoming quadrants to interleave with communication */
      p4est_balance_response (p4est, peer, btype, borders);
      qcount = peer->send_second.elem_count;
      if (qcount > 0) {
        P4EST_LDEBUGF ("Balance B send %llu quadrants to %d\n",
                       (unsigned long long) qcount, j);
        ++send_load[1];
      }
      else {
        ++send_zero[1];
      }
      peer->send_second_count = (int) qcount;
      if (qcount > 0) {
#ifdef P4EST_ENABLE_DEBUG
        checksum =
          p4est_quadrant_checksum (&peer->send_second, &checkarray, 0);
        P4EST_LDEBUGF ("Balance B send checksum 0x%08x to %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */
      }
      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_second.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_SECOND_LOAD,
                          p4est->mpicomm, &send_requests_second_load[j]);
      SC_CHECK_MPI (mpiret);
      ++request_second_send;
    }
  }
  for (j = 0; j < num_procs; ++j) {
    /* we cannot have exited this loop unless the first recv is complete */
    P4EST_ASSERT (recv_requests_first_load[j] == MPI_REQUEST_NULL);
  }
#endif /* P4EST_ENABLE_MPI */

  /* simulate send and receive with myself across tree boundaries */
  peer = peers + rank;
  sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
  qcount = peer->send_first.elem_count;
  peer->recv_first_count = peer->send_first_count = (int) qcount;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_first;
  sc_array_resize (qarray, qcount);
  memcpy (qarray->array, peer->send_first.array, qbytes);
  p4est_balance_response (p4est, peer, btype, borders);
  qcount = peer->send_second.elem_count;
  peer->recv_second_count = peer->send_second_count = (int) qcount;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_second;
  sc_array_resize (qarray, qcount);
  memcpy (qarray->array, peer->send_second.array, qbytes);

#ifdef P4EST_ENABLE_MPI
  /* receive second round appending to the same receive buffer */
  for (k = 0; k < request_first_send; k++) {
    MPI_Status          status;

    mpiret =
      MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_BALANCE_SECOND_LOAD,
                 p4est->mpicomm, &status);
    SC_CHECK_MPI (mpiret);
    j = status.MPI_SOURCE;
    peer = peers + j;
    mpiret = MPI_Get_count (&status, MPI_BYTE, &rcount);
    qcount = rcount / sizeof (p4est_quadrant_t);
    sc_array_resize (&(peer->recv_second), (size_t) qcount);
    total_recv_count += qcount;
    if (qcount) {
      recv_load[1]++;
    }
    else {
      recv_zero[1]++;
    }
    peer->recv_second_count = qcount;
    mpiret = MPI_Irecv (peer->recv_second.array, rcount, MPI_BYTE,
                        j, P4EST_COMM_BALANCE_SECOND_LOAD, p4est->mpicomm,
                        &recv_requests_second_load[j]);
    SC_CHECK_MPI (mpiret);
  }
  /* we cannot have initiated the second receive unless the first send is
   * complete */
  /* wait for the second recv to finish */
  mpiret = MPI_Waitall (num_procs, send_requests_first_load,
                        MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  /* wait for the second recv to finish */
  mpiret = MPI_Waitall (num_procs, recv_requests_second_load,
                        MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  for (j = 0; j < num_procs; ++j) {
    P4EST_ASSERT (recv_requests_second_load[j] == MPI_REQUEST_NULL);
#ifdef P4EST_ENABLE_DEBUG
    peer = peers + j;
    if (peer->recv_second.elem_count) {
      checksum = p4est_quadrant_checksum (&peer->recv_second, &checkarray, 0);
      P4EST_LDEBUGF ("Balance B recv checksum 0x%08x from %d\n", checksum, j);
    }
#endif
  }

  /* print buffer statistics */
  P4EST_VERBOSEF ("first send Z %d L %d recv Z %d L %d\n",
                  send_zero[0], send_load[0], recv_zero[0], recv_load[0]);
  P4EST_VERBOSEF ("second send Z %d L %d recv Z %d L %d\n",
                  send_zero[1], send_load[1], recv_zero[1], recv_load[1]);
  P4EST_VERBOSEF ("total send %d recv %d\n", total_send_count,
                  total_recv_count);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    if (peer->send_first_count > 0 || peer->recv_first_count > 0 ||
        peer->send_second_count > 0 || peer->recv_second_count > 0) {
      P4EST_VERBOSEF ("peer %d first S %llu R %d second S %llu R %d\n",
                      j, (unsigned long long) peer->send_first_count,
                      peer->recv_first_count,
                      (unsigned long long) peer->send_second_count,
                      peer->recv_second_count);
    }
  }
#endif /* P4EST_ENABLE_MPI */

  /* end balance_comm, start balance_B */
  if (p4est->inspect != NULL) {
    p4est->inspect->balance_comm += sc_MPI_Wtime ();
    p4est->inspect->balance_B = -sc_MPI_Wtime ();
    p4est->inspect->balance_B_count_in = 0;
    p4est->inspect->balance_B_count_out = 0;
    p4est->inspect->use_B = 1;
#ifdef P4EST_ENABLE_MPI
    for (k = 0; k < 2; ++k) {
      p4est->inspect->balance_load_sends[k] = send_load[k];
      p4est->inspect->balance_load_receives[k] = recv_load[k];
      p4est->inspect->balance_zero_sends[k] = send_zero[k];
      p4est->inspect->balance_zero_receives[k] = recv_zero[k];
    }
#endif
  }
  if (own_notify) {
    sc_notify_destroy (notify);
  }

  /* merge received quadrants */
  for (j = 0; j < num_procs; ++j) {
    size_t              fcount;

    /* access peer information */
    peer = peers + j;
    fcount = peer->recv_first.elem_count;
    qcount = fcount + peer->recv_second.elem_count;
    P4EST_ASSERT (peer->send_first_count ==
                  (int) peer->send_first.elem_count);
    P4EST_ASSERT (peer->send_second_count ==
                  (int) peer->send_second.elem_count);
    P4EST_ASSERT (peer->recv_second_count ==
                  (int) peer->recv_second.elem_count);
    if (qcount == 0) {
      continue;
    }

    /* merge received quadrants into correct tree */
    for (zz = 0; zz < qcount; ++zz) {
      s = zz < fcount ? p4est_quadrant_array_index (&peer->recv_first, zz) :
        p4est_quadrant_array_index (&peer->recv_second, zz - fcount);
      P4EST_ASSERT (p4est_quadrant_is_extended (s));
      qtree = s->p.piggy2.which_tree;
      if (qtree < first_tree || qtree > last_tree) {
        /* this is a corner/edge quadrant from the second pass of balance */
        continue;
      }
      if (borders == NULL) {
        tree = p4est_tree_array_index (p4est->trees, qtree);
        q = p4est_quadrant_array_push (&tree->quadrants);
        *q = *s;
        ++tree->quadrants_per_level[q->level];
        tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
        ++p4est->local_num_quadrants;
        p4est_quadrant_init_data (p4est, qtree, q, init_fn);
      }
      else {
        qarray = (sc_array_t *) sc_array_index (borders,
                                                (int) (qtree - first_tree));
        q = p4est_quadrant_array_push (qarray);
        *q = *s;
      }
    }
  }

  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    /* check if we are the only processor in an isolated tree */
    tree = p4est_tree_array_index (p4est->trees, nt);
    tree->quadrants_offset = p4est->local_num_quadrants;
    tquadrants = &tree->quadrants;
    treecount = tquadrants->elem_count;
    if (!(tree_flags[nt] & fully_owned_flag) ||
        (tree_flags[nt] & any_face_flag)) {
      /* we have most probably received quadrants, run sort and balance */
      /* balance the border, add it back into the tree, and linearize */
      p4est_balance_border (p4est, btype, nt, init_fn, replace_fn, borders);
      P4EST_VERBOSEF ("Balance tree %lld B %llu to %llu\n",
                      (long long) nt,
                      (unsigned long long) treecount,
                      (unsigned long long) tquadrants->elem_count);
    }
    p4est->local_num_quadrants += tquadrants->elem_count;
    tquadrants = NULL;          /* safeguard */
  }
  if (last_tree >= 0) {
    for (; nt < conn->num_trees; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  /* end balance_B */
  if (p4est->inspect != NULL) {
    p4est->inspect->balance_B += sc_MPI_Wtime ();
  }

#ifdef P4EST_ENABLE_MPI
  /* wait for all send operations */
  if (request_second_send > 0) {
    mpiret = MPI_Waitall (num_procs,
                          send_requests_second_load, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* compute global sum of send and receive counts */
#ifdef P4EST_ENABLE_DEBUG
  gtotal[0] = gtotal[1] = 0;
  ltotal[0] = (p4est_gloidx_t) total_send_count;
  ltotal[1] = (p4est_gloidx_t) total_recv_count;
  mpiret = MPI_Reduce (ltotal, gtotal, 2, P4EST_MPI_GLOIDX,
                       MPI_SUM, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_STATISTICSF ("Global number of shipped quadrants %lld\n",
                            (long long) gtotal[0]);
  P4EST_ASSERT (rank != 0 || gtotal[0] == gtotal[1]);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* loop over all local trees to finalize balance */
  all_outcount = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    all_outcount += tree->quadrants.elem_count;

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done balance tree %lld now %llu\n", (long long) nt,
                    (unsigned long long) tree->quadrants.elem_count);
  }

  /* cleanup temporary storage */
  P4EST_FREE (tree_flags);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    sc_array_reset (&peer->send_first);
    sc_array_reset (&peer->send_second);
    sc_array_reset (&peer->recv_first);
    sc_array_reset (&peer->recv_second);
  }
  P4EST_FREE (peers);

  if (borders != NULL) {
    for (zz = 0; zz < localcount; zz++) {
      qarray = (sc_array_t *) sc_array_index (borders, zz);
      sc_array_reset (qarray);
    }
    sc_array_destroy (borders);
  }

#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);

#ifdef P4EST_ENABLE_MPI
  sc_array_destroy (receivers);
  sc_array_destroy (senders);
  P4EST_FREE (send_requests_first_load);
  P4EST_FREE (recv_statuses);
  P4EST_FREE (wait_indices);
#ifdef P4EST_ENABLE_DEBUG
  sc_array_reset (&checkarray);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* some sanity checks */
  P4EST_ASSERT ((p4est_locidx_t) all_outcount == p4est->local_num_quadrants);
  P4EST_ASSERT (pre_adapt_flags || all_outcount >= all_incount);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + all_outcount - all_incount ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_VERBOSEF ("Balance skipped %lld\n", (long long) skipped);

  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_partition (p4est_t * p4est, int allow_for_coarsening,
                 p4est_weight_t weight_fn)
{
  (void) p4est_partition_ext (p4est, allow_for_coarsening, weight_fn);
}

/* p4est may not be valid yet: we expect
   each process to know its own local number of quadrants, but the global
   total may not be known yet */
p4est_locidx_t     *
p4est_partition_compute (p4est_t * p4est, int partition_for_coarsening,
                         p4est_gloidx_t global_num_quadrants,
                         p4est_weight_t weight_fn)
{
  p4est_gloidx_t      global_shipped = 0;
  sc_flopinfo_t       snap;
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
  int                 low_source, high_source;
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  const p4est_locidx_t local_num_quadrants = p4est->local_num_quadrants;
  int                 i, p;
  int                 send_lowest, send_highest;
  int                 num_sends, rcount, base_index;
  size_t              lz;
  ssize_t             lowers;
  p4est_topidx_t      nt;
  p4est_locidx_t      kl, qlocal;
  p4est_locidx_t     *num_quadrants_in_proc;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      send_index, recv_low, recv_high, qcount;
  p4est_gloidx_t     *send_array;
  int64_t             weight, weight_sum;
  int64_t             cut, my_lowcut, my_highcut;
  int64_t            *local_weights;    /* cumulative weights by quadrant */
  int64_t            *global_weight_sums;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  MPI_Request        *send_requests, recv_requests[2];
  MPI_Status          recv_statuses[2];
  p4est_gloidx_t      num_corrected;
#endif /* P4EST_ENABLE_MPI */

  P4EST_FUNC_SNAP (p4est, &snap);

  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc = P4EST_ALLOC (p4est_locidx_t, num_procs);

  if (weight_fn == NULL) {
    /* Divide up the quadrants equally */
    if (global_num_quadrants < 0) {
      /* global count is not known */
      p4est_gloidx_t      glocal_num_quadrants =
        (p4est_gloidx_t) local_num_quadrants;

      mpiret =
        sc_MPI_Allreduce (&glocal_num_quadrants, &global_num_quadrants, 1,
                          P4EST_MPI_GLOIDX, MPI_SUM, p4est->mpicomm);
      SC_CHECK_MPI (mpiret);
      P4EST_ASSERT (global_num_quadrants >= 0);
    }
    for (p = 0, next_quadrant = 0; p < num_procs; ++p) {
      prev_quadrant = next_quadrant;
      next_quadrant =
        p4est_partition_cut_gloidx (global_num_quadrants, p + 1, num_procs);
      qcount = next_quadrant - prev_quadrant;
      P4EST_ASSERT (0 <= qcount
                    && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
      num_quadrants_in_proc[p] = (p4est_locidx_t) (qcount);
    }
  }
  else {
    /* do a weighted partition */
    local_weights = P4EST_ALLOC (int64_t, local_num_quadrants + 1);
    global_weight_sums = P4EST_ALLOC (int64_t, num_procs + 1);
    P4EST_VERBOSEF ("local quadrant count %lld\n",
                    (long long) local_num_quadrants);

    /* linearly sum weights across all trees */
    kl = 0;
    local_weights[0] = 0;
    for (nt = first_tree; nt <= last_tree; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      for (lz = 0; lz < tree->quadrants.elem_count; ++lz, ++kl) {
        q = p4est_quadrant_array_index (&tree->quadrants, lz);
        weight = (int64_t) weight_fn (p4est, nt, q);
        P4EST_ASSERT (weight >= 0);
        local_weights[kl + 1] = local_weights[kl] + weight;
      }
    }
    P4EST_ASSERT (kl == local_num_quadrants);
    weight_sum = local_weights[local_num_quadrants];
    P4EST_VERBOSEF ("local weight sum %lld\n", (long long) weight_sum);

    /* distribute local weight sums */
    global_weight_sums[0] = 0;
    mpiret = MPI_Allgather (&weight_sum, 1, MPI_LONG_LONG_INT,
                            &global_weight_sums[1], 1, MPI_LONG_LONG_INT,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    /* adjust all arrays to reflect the global weight */
    for (i = 0; i < num_procs; ++i) {
      global_weight_sums[i + 1] += global_weight_sums[i];
    }
    if (rank > 0) {
      weight_sum = global_weight_sums[rank];
      for (kl = 0; kl <= local_num_quadrants; ++kl) {
        local_weights[kl] += weight_sum;
      }
    }
    P4EST_ASSERT (local_weights[0] == global_weight_sums[rank]);
    P4EST_ASSERT (local_weights[local_num_quadrants] ==
                  global_weight_sums[rank + 1]);
    weight_sum = global_weight_sums[num_procs];

    if (rank == 0) {
      for (i = 0; i <= num_procs; ++i) {
        P4EST_GLOBAL_VERBOSEF ("Global weight sum [%d] %lld\n",
                               i, (long long) global_weight_sums[i]);
      }
    }

    /* if all quadrants have zero weight we do nothing */
    if (weight_sum == 0) {
      P4EST_FREE (local_weights);
      P4EST_FREE (global_weight_sums);
      P4EST_FREE (num_quadrants_in_proc);

      /* in particular, there is no need to bumb the revision counter */
      P4EST_ASSERT (global_shipped == 0);

      P4EST_FUNC_SHOT (p4est, &snap);
      return NULL;              /* NULL inidicates no shipping to be done */
    }

    /* determine processor ids to send to */
    send_lowest = num_procs;
    send_highest = 0;
    for (i = 1; i <= num_procs; ++i) {
      cut = p4est_partition_cut_uint64 (weight_sum, i, num_procs);
      if (global_weight_sums[rank] < cut &&
          cut <= global_weight_sums[rank + 1]) {
        send_lowest = SC_MIN (send_lowest, i);
        send_highest = SC_MAX (send_highest, i);
      }
    }
    /*
     * send low cut to send_lowest..send_highest
     * and high cut to send_lowest-1..send_highest-1
     */
    P4EST_LDEBUGF ("my send peers %d %d\n", send_lowest, send_highest);

    num_sends = 2 * (send_highest - send_lowest + 1);
    if (num_sends <= 0) {
      num_sends = 0;
      send_requests = NULL;
      send_array = NULL;
    }
    else {
      send_requests = P4EST_ALLOC (MPI_Request, num_sends);
      send_array = P4EST_ALLOC (p4est_gloidx_t, num_sends);
      lowers = 0;
      for (i = send_lowest; i <= send_highest; ++i) {
        base_index = 2 * (i - send_lowest);
        if (i < num_procs) {
          /* do binary search in the weight array */
          lowers = sc_search_lower_bound64 (p4est_partition_cut_uint64
                                            (weight_sum, i, num_procs),
                                            local_weights,
                                            (size_t) local_num_quadrants + 1,
                                            (size_t) lowers);
          P4EST_ASSERT (lowers > 0
                        && (p4est_locidx_t) lowers <= local_num_quadrants);

          /* send low bound */
          send_index = send_array[base_index + 1] =
            (p4est_gloidx_t) lowers + p4est->global_first_quadrant[rank];
          P4EST_LDEBUGF ("send A %d %d %d index %lld base %d to %d\n",
                         send_lowest, i, send_highest,
                         (long long) send_index, base_index + 1, i);
          mpiret =
            MPI_Isend (&send_array[base_index + 1], 1, P4EST_MPI_GLOIDX, i,
                       P4EST_COMM_PARTITION_WEIGHTED_LOW, p4est->mpicomm,
                       &send_requests[base_index + 1]);
          SC_CHECK_MPI (mpiret);
        }
        else {
          lowers = 0;
          send_index = global_num_quadrants;
          send_requests[base_index + 1] = MPI_REQUEST_NULL;
          send_array[base_index + 1] = -1;
        }

        /* send high bound */
        send_array[base_index] = send_index;
        P4EST_LDEBUGF ("send B %d %d %d index %lld base %d to %d\n",
                       send_lowest, i, send_highest,
                       (long long) send_index, base_index, i - 1);
        mpiret = MPI_Isend (&send_array[base_index], 1, P4EST_MPI_GLOIDX,
                            i - 1, P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                            p4est->mpicomm, &send_requests[base_index]);
        SC_CHECK_MPI (mpiret);
      }
    }

    /* determine processor ids to receive from and post irecv */
    i = 0;
    my_lowcut = p4est_partition_cut_uint64 (weight_sum, rank, num_procs);
    if (my_lowcut == 0) {
      recv_low = 0;
      recv_requests[0] = MPI_REQUEST_NULL;
      low_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_lowcut &&
            my_lowcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("recv A from %d\n", i);
          mpiret = MPI_Irecv (&recv_low, 1, P4EST_MPI_GLOIDX, i,
                              P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm, &recv_requests[0]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      low_source = i;
    }
    my_highcut = p4est_partition_cut_uint64 (weight_sum, rank + 1, num_procs);
    if (my_highcut == 0) {
      recv_high = 0;
      recv_requests[1] = MPI_REQUEST_NULL;
      high_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_highcut &&
            my_highcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("recv B from %d\n", i);
          mpiret = MPI_Irecv (&recv_high, 1, P4EST_MPI_GLOIDX, i,
                              P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                              p4est->mpicomm, &recv_requests[1]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      high_source = i;
    }
    P4EST_LDEBUGF ("my recv peers %d %d cuts %lld %lld\n",
                   low_source, high_source,
                   (long long) my_lowcut, (long long) my_highcut);

    /* free temporary memory */
    P4EST_FREE (local_weights);
    P4EST_FREE (global_weight_sums);

    /* wait for sends and receives to complete */
    if (num_sends > 0) {
      mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
      SC_CHECK_MPI (mpiret);
      P4EST_FREE (send_requests);
      P4EST_FREE (send_array);
    }
    mpiret = MPI_Waitall (2, recv_requests, recv_statuses);
    SC_CHECK_MPI (mpiret);
    if (my_lowcut != 0) {
      SC_CHECK_ABORT (recv_statuses[0].MPI_SOURCE == low_source,
                      "Wait low source");
      SC_CHECK_ABORT (recv_statuses[0].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_LOW, "Wait low tag");
      mpiret = MPI_Get_count (&recv_statuses[0], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait low count %d", rcount);
    }
    if (my_highcut != 0) {
      SC_CHECK_ABORT (recv_statuses[1].MPI_SOURCE == high_source,
                      "Wait high source");
      SC_CHECK_ABORT (recv_statuses[1].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_HIGH, "Wait high tag");
      mpiret = MPI_Get_count (&recv_statuses[1], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait high count %d", rcount);
    }

    /* communicate the quadrant ranges */
    qcount = recv_high - recv_low;
    P4EST_LDEBUGF ("weighted partition count %lld\n", (long long) qcount);
    P4EST_ASSERT (qcount >= 0 && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
    qlocal = (p4est_locidx_t) qcount;
    mpiret = MPI_Allgather (&qlocal, 1, P4EST_MPI_LOCIDX,
                            num_quadrants_in_proc, 1, P4EST_MPI_LOCIDX,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

#if(0)
    /* run through the count array and repair zero ranges */
    for (i = 0; i < num_procs; ++i) {
      if (num_quadrants_in_proc[i] == 0) {
        for (p = i - 1; p >= 0; --p) {
          P4EST_ASSERT (num_quadrants_in_proc[p] > 0);
          if (num_quadrants_in_proc[p] > 1) {
            --num_quadrants_in_proc[p];
            ++num_quadrants_in_proc[i];
            break;
          }
        }
        if (p < 0) {
          for (p = i + 1; p < num_procs; ++p) {
            P4EST_ASSERT (num_quadrants_in_proc[p] >= 0);
            if (num_quadrants_in_proc[p] > 1) {
              --num_quadrants_in_proc[p];
              ++num_quadrants_in_proc[i];
              break;
            }
          }
          P4EST_ASSERT (p < num_procs);
        }
      }
    }
#endif
  }

  /* correct partition */
  if (partition_for_coarsening) {
    num_corrected =
      p4est_partition_for_coarsening (p4est, num_quadrants_in_proc);
    P4EST_GLOBAL_INFOF
      ("Designated partition for coarsening %lld quadrants moved\n",
       (long long) num_corrected);
  }

  P4EST_FUNC_SHOT (p4est, &snap);
  return num_quadrants_in_proc;
}

p4est_gloidx_t
p4est_partition_ext (p4est_t * p4est, int partition_for_coarsening,
                     p4est_weight_t weight_fn)
{
  p4est_gloidx_t      global_shipped = 0;
  const p4est_gloidx_t global_num_quadrants = p4est->global_num_quadrants;
#ifdef P4EST_ENABLE_MPI
  p4est_locidx_t     *num_quadrants_in_proc;
#endif /* P4EST_ENABLE_MPI */
  sc_flopinfo_t       snap;

  P4EST_FUNC_SNAP (p4est, &snap);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF
    ("Into " P4EST_STRING
     "_partition with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

  /* this function does nothing in a serial setup */
  if (p4est->mpisize == 1) {
    P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_partition no shipping\n");

    /* in particular, there is no need to bumb the revision counter */
    P4EST_ASSERT (global_shipped == 0);
    return global_shipped;
  }

  p4est_log_indent_push ();

#ifdef P4EST_ENABLE_MPI
  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc =
    p4est_partition_compute (p4est, partition_for_coarsening,
                             p4est->global_num_quadrants, weight_fn);

  if (!num_quadrants_in_proc) {
    P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_partition no shipping\n");
    P4EST_FUNC_SHOT (p4est, &snap);
    return 0;
  }

  /* run the partition algorithm with proper quadrant counts */
  global_shipped = p4est_partition_given (p4est, num_quadrants_in_proc);
  if (global_shipped) {
    /* the partition of the forest has changed somewhere */
    ++p4est->revision;
  }
  P4EST_FREE (num_quadrants_in_proc);

  /* check validity of the p4est */
  P4EST_ASSERT (p4est_is_valid (p4est));
#endif /* P4EST_ENABLE_MPI */

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_partition shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / global_num_quadrants);

  P4EST_FUNC_SHOT (p4est, &snap);
  return global_shipped;
}

p4est_gloidx_t
p4est_partition_for_coarsening (p4est_t * p4est,
                                p4est_locidx_t * num_quadrants_in_proc)
{
#ifdef P4EST_ENABLE_MPI
  int                 num_procs = p4est->mpisize;
  int                 rank = p4est->mpirank;
  int                 mpiret;
  p4est_gloidx_t      global_num_quadrants = p4est->global_num_quadrants;
  int                 i, send_lowest, send_highest, num_sends;
  int                 parent_index;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  p4est_locidx_t      num_quadrants_in_tree;
  p4est_topidx_t      it, tree_index;
  p4est_gloidx_t      iq, quad_id_near_cut;
  p4est_gloidx_t      min_quad_id, max_quad_id;
  int8_t              quad_near_cut_level;
  p4est_gloidx_t     *partition_now = p4est->global_first_quadrant;
  p4est_gloidx_t     *partition_new;
  p4est_quadrant_t   *parent_send;
  MPI_Request        *send_requests;
  MPI_Request        *receive_requests;
  int                 receive_lowest, receive_highest, num_receives;
  int                 process_with_cut = -1, process_with_cut_recv_id = -1;
  p4est_quadrant_t   *parent_receive;
  int                *receive_process;
  int                *correction, correction_local = 0;
  int                 current_proc, next_proc;
  p4est_gloidx_t      num_moved_quadrants;
  sc_flopinfo_t       snap;

  P4EST_FUNC_SNAP (p4est, &snap);

  /* create array with first quadrants of new partition */
  partition_new = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  partition_new[0] = 0;
  for (i = 1; i < num_procs; i++) {     /* loop over all processes */
    partition_new[i] = partition_new[i - 1] + num_quadrants_in_proc[i - 1];
  }
  partition_new[num_procs] = global_num_quadrants;

  /* BEGIN: send */
  if (partition_now[rank] < partition_now[rank + 1]) {
    /* if this process has quadrants */
    /* determine number and min/max process ids to send to */
    num_sends = 0;              /* number of sends */
    send_lowest = num_procs;    /* lowest process id */
    send_highest = 0;           /* highest process id */
    for (i = 1; i < num_procs; i++) {
      /* loop over all processes (without first) */
      if (partition_new[i] < partition_new[i + 1] &&
          partition_now[rank] <= partition_new[i] + P4EST_CHILDREN - 2 &&
          partition_new[i] - P4EST_CHILDREN + 1 < partition_now[rank + 1]) {
        /* if this process has relevant quadrants for process `i` */
        num_sends++;
        send_lowest = SC_MIN (send_lowest, i);
        send_highest = SC_MAX (send_highest, i);
      }
    }
  }
  else {
    /* set number of messages to send */
    num_sends = 0;
  }

  if (num_sends > 0) {          /* if this process sends messages */
    /* allocate send messages */
    send_requests = P4EST_ALLOC (MPI_Request, num_sends);
    parent_send = P4EST_ALLOC (p4est_quadrant_t, num_sends);

    /* array index of send messages */
    parent_index = 0;

    for (i = send_lowest; i <= send_highest; i++) {
      /* loop over all process candidates to send to */
      if (!(partition_new[i] < partition_new[i + 1] &&
            partition_now[rank] <= partition_new[i] + P4EST_CHILDREN - 2 &&
            partition_new[i] - P4EST_CHILDREN + 1 <
            partition_now[rank + 1])) {
        /* if this process has no relevant quadrants for process `i` */
        continue;
      }

      /* get nearest quadrant `quad_id_near_cut` to cut `partition_new[i]` */
      if (partition_now[rank] <= partition_new[i] &&
          partition_new[i] < partition_now[rank + 1]) {
        /* if cut is owned by this process */
        quad_id_near_cut = partition_new[i];
      }
      else {
        if (P4EST_GLOIDX_ABS (partition_new[i] - partition_now[rank]) <
            P4EST_GLOIDX_ABS (partition_new[i] - partition_now[rank + 1] +
                              1)) {
          quad_id_near_cut = partition_now[rank];
        }
        else {
          quad_id_near_cut = partition_now[rank + 1] - 1;
        }
      }

      /* get tree `tree` of quadrant `quad_id_near_cut` */
      num_quadrants_in_tree = partition_now[rank + 1] - partition_now[rank];
      tree_index = -1;
      for (it = p4est->first_local_tree; it <= p4est->last_local_tree; it++) {
        /* loop over all local trees */
        tree = p4est_tree_array_index (p4est->trees, it);
        if (tree->quadrants_offset <= quad_id_near_cut - partition_now[rank]) {
          tree_index = it;
        }
        else {
          num_quadrants_in_tree = tree->quadrants_offset;
          break;
        }
      }
      tree = p4est_tree_array_index (p4est->trees, tree_index);
      num_quadrants_in_tree -= tree->quadrants_offset;

      /* get quadrant with index `quad_id_near_cut` */
      q = p4est_quadrant_array_index (&tree->quadrants,
                                      quad_id_near_cut - partition_now[rank] -
                                      tree->quadrants_offset);

      /* get level of quadrant near cut */
      quad_near_cut_level = q->level;

      if (quad_near_cut_level > 0) {
        /* if quadrant near cut is not root of tree, i.e. level is not zero */
        /* get parent of quadrant near cut */
        p4est_quadrant_parent (q, &parent_send[parent_index]);

        /* get min quadrant with same parent */
        min_quad_id = quad_id_near_cut;
        for (iq = quad_id_near_cut;
             iq >= SC_MAX (partition_now[rank] + tree->quadrants_offset,
                           partition_new[i] - P4EST_CHILDREN + 1); iq--) {
          /* loop over eligible quadrants */
          /* get quadrant with index `iq` */
          q = p4est_quadrant_array_index (&tree->quadrants,
                                          iq - partition_now[rank] -
                                          tree->quadrants_offset);

          /* check quadrant `iq` */
          if (q->level == quad_near_cut_level) {        /* if same level */
            if (p4est_quadrant_is_parent (&parent_send[parent_index], q)) {
              /* if same parent */
              min_quad_id = iq;
            }
            else {
              break;
            }
          }
          else {
            break;
          }
        }

        /* get max quadrant with same parent */
        max_quad_id = quad_id_near_cut;
        for (iq = quad_id_near_cut;
             iq <=
             SC_MIN (partition_now[rank] + tree->quadrants_offset +
                     num_quadrants_in_tree - 1,
                     partition_new[i] + P4EST_CHILDREN - 2); iq++) {
          /* loop over eligible quadrants */

          /* get quadrant `iq` */
          q = p4est_quadrant_array_index (&tree->quadrants,
                                          iq - partition_now[rank] -
                                          tree->quadrants_offset);

          /* check quadrant `iq` */
          if (q->level == quad_near_cut_level) {        /* if same level */
            if (p4est_quadrant_is_parent (&parent_send[parent_index], q)) {
              /* if same parent */
              max_quad_id = iq;
            }
            else {
              break;
            }
          }
          else {
            break;
          }
        }

        /* write tree */
        parent_send[parent_index].p.piggy3.which_tree = tree_index;

        if (quad_id_near_cut == partition_new[i]) {
          /* if this process has cut */
          /* encode number of quadrants with same parent before and after
           * `partition_new[i]` into one integer */
          parent_send[parent_index].p.piggy3.local_num =
            (partition_new[i] - min_quad_id) * P4EST_CHILDREN +
            (max_quad_id - partition_new[i]);
        }
        else {
          /* write number of quadrants with same parent */
          parent_send[parent_index].p.piggy3.local_num =
            max_quad_id - min_quad_id + 1;
        }

        /* MPI send: parent */
        mpiret = MPI_Isend (&parent_send[parent_index],
                            sizeof (p4est_quadrant_t), MPI_BYTE,
                            i,
                            P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                            &send_requests[parent_index]);
        SC_CHECK_MPI (mpiret);
      }
      else {
        /* if quadrant near cut is root of tree, i.e., level is zero,
           set parent as tree root `q` */
        parent_send[parent_index].level = q->level;
        parent_send[parent_index].x = q->x;
        parent_send[parent_index].y = q->y;
#ifdef P4_TO_P8
        parent_send[parent_index].z = q->z;
#endif

        /* write tree */
        parent_send[parent_index].p.piggy3.which_tree = tree_index;

        /* write number of quadrants with same "parent" */
        parent_send[parent_index].p.piggy3.local_num = 1;

        /* MPI send: root of tree */
        mpiret = MPI_Isend (&parent_send[parent_index],
                            sizeof (p4est_quadrant_t), MPI_BYTE,
                            i,
                            P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                            &send_requests[parent_index]);
        SC_CHECK_MPI (mpiret);
      }

      /* increment parent index */
      parent_index++;
    }
  }
  /* END: send */

  /* BEGIN: receive */
  if (rank != 0 && partition_new[rank] < partition_new[rank + 1]) {
    /* if this process should get quadrants */
    /* determine process ids to receive from */
    num_receives = 0;           /* number of receives */
    receive_lowest = num_procs; /* lowest process id */
    receive_highest = 0;        /* highest process id */
    for (i = 0; i < num_procs; i++) {   /* loop over all processes */
      if (partition_now[i] < partition_now[i + 1] &&
          partition_now[i] <= partition_new[rank] + P4EST_CHILDREN - 2 &&
          partition_new[rank] - P4EST_CHILDREN + 1 < partition_now[i + 1]) {
        /* if process `i` has relevant quadrants for this process */
        num_receives++;
        receive_lowest = SC_MIN (receive_lowest, i);
        receive_highest = SC_MAX (receive_highest, i);

        if (partition_now[i] <= partition_new[rank] &&
            partition_new[rank] < partition_now[i + 1]) {
          /* if cut is owned by process `i` */
          /* process id that sends parent of cut quadrant */
          process_with_cut = i;

          /* array index of receive messages of process with cut  */
          process_with_cut_recv_id = num_receives - 1;
        }
      }
    }
  }
  else {
    /* set number of messages to receive */
    num_receives = 0;

    /* set correction */
    correction_local = 0;
  }

  if (num_receives > 0) {       /* if this process receives messages */
    /* allocate receive messages */
    receive_requests = P4EST_ALLOC (MPI_Request, num_receives);
    parent_receive = P4EST_ALLOC (p4est_quadrant_t, num_receives);
    receive_process = P4EST_ALLOC (int, num_receives);

    /* array index of receive messages */
    parent_index = 0;

    for (i = receive_lowest; i <= receive_highest; i++) {
      /* loop over all process candidates to receive from */
      if (!(partition_now[i] < partition_now[i + 1] &&
            partition_now[i] <= partition_new[rank] + P4EST_CHILDREN - 2 &&
            partition_new[rank] - P4EST_CHILDREN + 1 <
            partition_now[i + 1])) {
        /* if process `i` has no relevant quadrants for this process */
        continue;
      }
      /* store process index */
      receive_process[parent_index] = i;

      /* MPI receive */
      mpiret = MPI_Irecv (&parent_receive[parent_index],
                          sizeof (p4est_quadrant_t), MPI_BYTE,
                          i,
                          P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                          &receive_requests[parent_index]);
      SC_CHECK_MPI (mpiret);

      /* increment parent index */
      parent_index++;
    }
  }
  /* END: receive */

  /* BEGIN: wait for MPI receive to complete */
  if (num_receives > 0) {
    /* wait for receives to complete */
    mpiret =
      MPI_Waitall (num_receives, receive_requests, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* free receive memory */
    P4EST_FREE (receive_requests);
  }
  /* END: wait for MPI recieve to complete */

  /* BEGIN: compute correction with received quadrants */
  if (num_receives > 0) {
    /* if this process received quadrants */
    min_quad_id = partition_new[rank];  /* min quadrant id with same parent */
    max_quad_id = partition_new[rank];  /* max quadrant id with same parent */

    for (i = 0; i < num_receives; i++) {
      /* loop over all received (parent or tree root) quadrants */
      if (parent_receive[i].p.piggy3.which_tree ==
          parent_receive[process_with_cut_recv_id].p.piggy3.which_tree
          &&
          p4est_quadrant_is_equal (&parent_receive[i],
                                   &parent_receive[process_with_cut_recv_id]
          )) {
        /* if trees and parents are equal */
        /* decrease/increase min/max quadrant with same parent */
        if (receive_process[i] < process_with_cut) {
          /* if before process with cut */
          /* decrease min quadrant */
          min_quad_id -= parent_receive[i].p.piggy3.local_num;
        }
        else if (receive_process[i] > process_with_cut) {
          /* if after process with cut */
          /* increase max quadrant */
          max_quad_id += parent_receive[i].p.piggy3.local_num;
        }
        else {
          /* decrease min quadrant */
          min_quad_id -=
            parent_receive[i].p.piggy3.local_num / P4EST_CHILDREN;

          /* increase max quadrant */
          max_quad_id +=
            parent_receive[i].p.piggy3.local_num % P4EST_CHILDREN;
        }
      }
    }

    /* compute correction */
    correction_local =
      (int) p4est_partition_correction (partition_new, num_procs, rank,
                                        min_quad_id, max_quad_id);

    /* free receive memory */
    P4EST_FREE (parent_receive);
    P4EST_FREE (receive_process);
  }
  /* END: compute correction with received parent quadrants */

  /* free memory */
  P4EST_FREE (partition_new);

  /* communicate corrections */
  correction = P4EST_ALLOC (int, num_procs);
  mpiret = MPI_Allgather (&correction_local, 1, MPI_INT,
                          correction, 1, MPI_INT, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* BEGIN: wait for MPI send to complete */
  if (num_sends > 0) {
    /* wait for sends to complete */
    mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* free send memory */
    P4EST_FREE (parent_send);
    P4EST_FREE (send_requests);
  }
  /* END: wait for MPI send to complete */

  /* correct partition */
  current_proc =
    p4est_next_nonempty_process (0, num_procs, num_quadrants_in_proc);
  next_proc =
    p4est_next_nonempty_process (current_proc + 1, num_procs,
                                 num_quadrants_in_proc);
  num_moved_quadrants = 0;
  while (current_proc < num_procs) {
    /* loop over all non empty processes */
    /* compute correct partition for process `current_proc` */
    if (0 < current_proc && current_proc < num_procs) {
      /* if any process but first */
      num_quadrants_in_proc[current_proc] += correction[current_proc];
      /* input is just a locidx, but the result is gloidx so we cast cleanly */
      num_moved_quadrants += P4EST_GLOIDX_ABS (correction[current_proc]);
    }
    if (next_proc < num_procs) {
      /* if first process or next process is feasible */
      num_quadrants_in_proc[current_proc] -= correction[next_proc];
    }

    /* increase process ids */
    current_proc = next_proc;
    next_proc =
      p4est_next_nonempty_process (next_proc + 1, num_procs,
                                   num_quadrants_in_proc);
  }

  /* free memory */
  P4EST_FREE (correction);

  /* return absolute number of moved quadrants */
  P4EST_FUNC_SHOT (p4est, &snap);
  return num_moved_quadrants;
#else
  P4EST_ASSERT (num_quadrants_in_proc[0] == p4est->local_num_quadrants);

  return 0;
#endif
}

unsigned
p4est_checksum (p4est_t * p4est)
{
#ifdef P4EST_HAVE_ZLIB
  uLong               treecrc, crc;
  size_t              scount, ssum;
  p4est_topidx_t      nt;
  p4est_tree_t       *tree;
  sc_array_t          checkarray;

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&checkarray, 4);
  crc = adler32 (0, Z_NULL, 0);
  ssum = 0;
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    treecrc =
      (uLong) p4est_quadrant_checksum (&tree->quadrants, &checkarray, 0);
    scount = 4 * checkarray.elem_count;
    ssum += scount;
    crc = adler32_combine (crc, treecrc, (z_off_t) scount);
  }
  sc_array_reset (&checkarray);
  P4EST_ASSERT ((p4est_locidx_t) ssum ==
                p4est->local_num_quadrants * 4 * (P4EST_DIM + 1));

  return p4est_comm_checksum (p4est, (unsigned) crc, ssum);
#else
  sc_abort_collective
    ("Configure did not find a recent enough zlib.  Abort.\n");

  return 0;
#endif /* !P4EST_HAVE_ZLIB */
}

void
p4est_save (const char *filename, p4est_t * p4est, int save_data)
{
  p4est_save_ext (filename, p4est, save_data, 1);
}

void
p4est_save_ext (const char *filename, p4est_t * p4est,
                int save_data, int save_partition)
{
  const int           headc = 6;
  const int           align = 32;
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
#ifndef P4EST_MPIIO_WRITE
  MPI_Status          mpistatus;
#endif
#endif
  int                 retval;
  int                 num_procs, save_num_procs, rank;
  int                 i;
  long                fpos = -1, foffset;
  size_t              data_size, qbuf_size, comb_size, head_count;
  size_t              zz, zcount;
  uint64_t           *u64a;
  FILE               *file;
#ifdef P4EST_MPIIO_WRITE
  MPI_File            mpifile;
  MPI_Offset          mpipos;
  MPI_Offset          mpithis;
#else
  long                fthis;
#endif
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t     *pertree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  char               *lbuf, *bp;
  p4est_qcoord_t     *qpos;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_save %s\n", filename);
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (p4est->connectivity));
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* when data is not saved the size is set to zero */
  data_size = save_data ? p4est->data_size : 0;

  /* zero data size is effectively not saved */
  if (data_size == 0) {
    save_data = 0;
  }

  /* other parameters */
  num_trees = p4est->connectivity->num_trees;
  num_procs = p4est->mpisize;
  save_num_procs = save_partition ? num_procs : 1;
  head_count = (size_t) (headc + save_num_procs) + (size_t) num_trees;
  rank = p4est->mpirank;
  qbuf_size = (P4EST_DIM + 1) * sizeof (p4est_qcoord_t);
  comb_size = qbuf_size + data_size;
  pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  p4est_comm_count_pertree (p4est, pertree);

  if (rank == 0) {
    p4est_connectivity_save (filename, p4est->connectivity);

    /* open file after writing connectivity to it */
    file = fopen (filename, "ab");
    SC_CHECK_ABORT (file != NULL, "file open");

    /* align the start of the header */
    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "first file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "first file align");
      ++fpos;
    }

    /* write format and partition information */
    u64a = P4EST_ALLOC (uint64_t, head_count);
    u64a[0] = P4EST_ONDISK_FORMAT;
    u64a[1] = (uint64_t) sizeof (p4est_qcoord_t);
    u64a[2] = (uint64_t) sizeof (p4est_quadrant_t);
    u64a[3] = (uint64_t) data_size;
    u64a[4] = (uint64_t) save_data;
    u64a[5] = (uint64_t) save_num_procs;
    if (save_partition) {
      P4EST_ASSERT (save_num_procs == num_procs);
      for (i = 0; i < num_procs; ++i) {
        u64a[headc + i] = (uint64_t) p4est->global_first_quadrant[i + 1];
      }
    }
    else {
      P4EST_ASSERT (save_num_procs == 1);
      u64a[headc] = (uint64_t) p4est->global_first_quadrant[num_procs];
    }
    for (jt = 0; jt < num_trees; ++jt) {
      u64a[headc + save_num_procs + jt] = (uint64_t) pertree[jt + 1];
    }
    sc_fwrite (u64a, sizeof (uint64_t), head_count,
               file, "write header information");
    P4EST_FREE (u64a);
    fpos += head_count * sizeof (uint64_t);

    /* align the start of the quadrants */
    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "second file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "second file align");
      ++fpos;
    }

#ifdef P4EST_MPIIO_WRITE
    /* We will close the sequential access to the file */
    /* best attempt to flush file to disk */
    retval = fflush (file);
    SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
    retval = fsync (fileno (file));
    SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
    retval = fclose (file);
    SC_CHECK_ABORT (retval == 0, "file close");
    file = NULL;
#else
    /* file is still open for sequential write mode */
#endif
  }
  else {
    file = NULL;
  }
  P4EST_FREE (pertree);

#ifndef P4EST_MPIIO_WRITE
  if (rank > 0) {
    /* wait for sequential synchronization */
#ifdef P4EST_ENABLE_MPI
    mpiret = MPI_Recv (&fpos, 1, MPI_LONG, rank - 1, P4EST_COMM_SAVE,
                       p4est->mpicomm, &mpistatus);
    SC_CHECK_MPI (mpiret);
#endif

    /* open file after all previous processors have written to it */
    file = fopen (filename, "rb+");
    SC_CHECK_ABORT (file != NULL, "file open");
  }
#else
  /* Every core opens the file in append mode */
  mpiret = MPI_File_open (p4est->mpicomm, (char *) filename,
                          MPI_MODE_WRONLY | MPI_MODE_APPEND |
                          MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &mpifile);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_File_get_position (mpifile, &mpipos);
  SC_CHECK_MPI (mpiret);
#endif

  if (rank > 0) {
    /* seek to the beginning of this processor's storage */
    foffset = (long) (p4est->global_first_quadrant[rank] * comb_size);

#ifndef P4EST_MPIIO_WRITE
    fthis = fpos + foffset;
    retval = fseek (file, fthis, SEEK_SET);
    SC_CHECK_ABORT (retval == 0, "seek data");
#else
    mpithis = mpipos + (MPI_Offset) foffset;
    mpiret = MPI_File_seek (mpifile, mpithis, MPI_SEEK_SET);
    SC_CHECK_MPI (mpiret);
#endif
  }

  /* write quadrant coordinates and data interleaved */
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
    zcount = tquadrants->elem_count;

    /* storage that will be written for this tree */
    bp = lbuf = P4EST_ALLOC (char, comb_size * zcount);
    for (zz = 0; zz < zcount; ++zz) {
      qpos = (p4est_locidx_t *) bp;
      q = p4est_quadrant_array_index (tquadrants, zz);
      *qpos++ = q->x;
      *qpos++ = q->y;
#ifdef P4_TO_P8
      *qpos++ = q->z;
#endif
      *qpos++ = (p4est_qcoord_t) q->level;
      if (save_data) {
        memcpy (qpos, q->p.user_data, data_size);
      }
      bp += comb_size;
    }
#ifndef P4EST_MPIIO_WRITE
    sc_fwrite (lbuf, comb_size, zcount, file, "write quadrants");
#else
    sc_mpi_write (mpifile, lbuf, comb_size * zcount, MPI_BYTE,
                  "write quadrants");
#endif
    P4EST_FREE (lbuf);
  }

#ifndef P4EST_MPIIO_WRITE
  /* best attempt to flush file to disk */
  retval = fflush (file);
  SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
  retval = fsync (fileno (file));
  SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
  retval = fclose (file);
  SC_CHECK_ABORT (retval == 0, "file close");
  file = NULL;

  /* initiate sequential synchronization */
#ifdef P4EST_ENABLE_MPI
  if (rank < num_procs - 1) {
    mpiret = MPI_Send (&fpos, 1, MPI_LONG, rank + 1, P4EST_COMM_SAVE,
                       p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
#endif
#else
  mpiret = MPI_File_close (&mpifile);
  SC_CHECK_MPI (mpiret);
#endif

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_save\n");
}

p4est_t            *
p4est_load (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
            int load_data, void *user_pointer,
            p4est_connectivity_t ** connectivity)
{
  return p4est_load_ext (filename, mpicomm, data_size, load_data,
                         0, 0, user_pointer, connectivity);
}

p4est_t            *
p4est_load_ext (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
                int load_data, int autopartition, int broadcasthead,
                void *user_pointer, p4est_connectivity_t ** connectivity)
{
  int                 retval;
  p4est_t            *p4est;
  sc_io_source_t     *src;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_load %s\n", filename);
  p4est_log_indent_push ();

  /* open file on all processors */

  src = sc_io_source_new (SC_IO_TYPE_FILENAME, SC_IO_ENCODE_NONE, filename);
  SC_CHECK_ABORT (src != NULL, "file source: possibly file not found");

  p4est = p4est_source_ext (src, mpicomm, data_size, load_data, autopartition,
                            broadcasthead, user_pointer, connectivity);

  retval = sc_io_source_destroy (src);
  SC_CHECK_ABORT (!retval, "source destroy");

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_load with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

  return p4est;
}

p4est_t            *
p4est_source_ext (sc_io_source_t * src, sc_MPI_Comm mpicomm, size_t data_size,
                  int load_data, int autopartition, int broadcasthead,
                  void *user_pointer, p4est_connectivity_t ** connectivity)
{
  const int           headc = 6;
  const int           align = 32;
  int                 root = 0;
  int                 retval;
  int                 mpiret;
  int                 num_procs, rank;
  int                 save_num_procs;
  int                 save_data;
  int                 i;
  uint64_t           *u64a, u64int;
  size_t              conn_bytes, file_offset;
  size_t              save_data_size;
  size_t              qbuf_size, comb_size, head_count;
  size_t              zz, zcount, zpadding;
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t     *gfq;
  p4est_gloidx_t     *pertree;
  p4est_qcoord_t     *qap;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  sc_array_t         *qarr, *darr;
  char               *dap, *lbuf;

  /* set some parameters */
  P4EST_ASSERT (src->bytes_out == 0);
  P4EST_ASSERT (connectivity != NULL);
  if (data_size == 0) {
    load_data = 0;
  }
  qbuf_size = (P4EST_DIM + 1) * sizeof (p4est_qcoord_t);

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* the first part of the header determines further offsets */
  save_data_size = (size_t) ULONG_MAX;
  save_num_procs = -1;
  conn = NULL;
  conn_bytes = 0;
  u64a = P4EST_ALLOC (uint64_t, headc + 1);
  if (!broadcasthead || rank == root) {

    /* read the forest connectivity */
    conn = p4est_connectivity_source (src);
    SC_CHECK_ABORT (conn != NULL, "connectivity source");
    zcount = src->bytes_out;
    zpadding = (align - zcount % align) % align;
    retval = sc_io_source_read (src, NULL, zpadding, NULL);
    SC_CHECK_ABORT (!retval, "source padding");
    conn_bytes = src->bytes_out;

    /* read format and some basic partition parameters */
    retval = sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t) headc,
                                NULL);
    SC_CHECK_ABORT (!retval, "read format");
    SC_CHECK_ABORT (u64a[0] == P4EST_ONDISK_FORMAT, "invalid format");
    SC_CHECK_ABORT (u64a[1] == (uint64_t) sizeof (p4est_qcoord_t),
                    "invalid coordinate size");
    SC_CHECK_ABORT (u64a[2] == (uint64_t) sizeof (p4est_quadrant_t),
                    "invalid quadrant size");
    save_data_size = (size_t) u64a[3];
    save_data = (int) u64a[4];
    if (load_data) {
      SC_CHECK_ABORT (save_data_size == data_size, "invalid data size");
      SC_CHECK_ABORT (save_data, "quadrant data not saved");
    }
    save_num_procs = (int) u64a[5];
    SC_CHECK_ABORT (autopartition || num_procs == save_num_procs,
                    "num procs mismatch");

    /* piggy-back the bytes for the connectivity onto first message */
    u64a[headc + 0] = (uint64_t) conn_bytes;
  }
  if (broadcasthead) {

    /* broadcast connectivity and first part of header */
    conn = p4est_connectivity_bcast (conn, root, mpicomm);
    mpiret = sc_MPI_Bcast (u64a, headc + 1, sc_MPI_LONG_LONG_INT,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
    if (rank != root) {

      /* make sure the rest of the processes has the information */
      SC_CHECK_ABORT (u64a[0] == P4EST_ONDISK_FORMAT, "invalid format");
      save_data_size = (size_t) u64a[3];
      save_data = (int) u64a[4];
      save_num_procs = (int) u64a[5];
      conn_bytes = (size_t) u64a[headc + 0];
    }
  }
  P4EST_ASSERT (save_num_procs >= 0);
  P4EST_ASSERT (save_data_size != (size_t) ULONG_MAX);
  *connectivity = conn;
  comb_size = qbuf_size + save_data_size;
  file_offset = conn_bytes + headc * sizeof (uint64_t);

  /* create partition data */
  gfq = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  gfq[0] = 0;
  if (!broadcasthead || rank == root) {
    if (!autopartition) {
      P4EST_ASSERT (num_procs == save_num_procs);
      u64a = P4EST_REALLOC (u64a, uint64_t, num_procs);
      sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t) num_procs,
                         NULL);
      SC_CHECK_ABORT (!retval, "read quadrant partition");
      for (i = 0; i < num_procs; ++i) {
        gfq[i + 1] = (p4est_gloidx_t) u64a[i];
      }
    }
    else {
      /* ignore saved partition and compute a new uniform one */
      retval = sc_io_source_read
        (src, NULL, (long) ((save_num_procs - 1) * sizeof (uint64_t)), NULL);
      SC_CHECK_ABORT (!retval, "seek over ignored partition");
      retval = sc_io_source_read (src, &u64int, sizeof (uint64_t), NULL);
      SC_CHECK_ABORT (!retval, "read quadrant count");
      for (i = 1; i <= num_procs; ++i) {
        gfq[i] = p4est_partition_cut_uint64 (u64int, i, num_procs);
      }
    }
  }
  if (broadcasthead) {
    mpiret = sc_MPI_Bcast (gfq + 1, num_procs, P4EST_MPI_GLOIDX,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  zcount = (size_t) (gfq[rank + 1] - gfq[rank]);
  file_offset += save_num_procs * sizeof (uint64_t);

  /* read pertree data */
  num_trees = conn->num_trees;
  pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  pertree[0] = 0;
  if (!broadcasthead || rank == root) {
    u64a = P4EST_REALLOC (u64a, uint64_t, num_trees);
    retval = sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t)
                                num_trees, NULL);
    SC_CHECK_ABORT (!retval, "read pertree information");
    for (jt = 0; jt < num_trees; ++jt) {
      pertree[jt + 1] = (p4est_gloidx_t) u64a[jt];
    }
    SC_CHECK_ABORT (gfq[num_procs] == pertree[num_trees], "pertree mismatch");
  }
  if (broadcasthead) {
    mpiret = sc_MPI_Bcast (pertree + 1, num_trees, P4EST_MPI_GLOIDX,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_FREE (u64a);
  file_offset += num_trees * sizeof (uint64_t);

  /* seek to the beginning of this processor's storage */
  if (!broadcasthead || rank == root) {
    P4EST_ASSERT (file_offset == src->bytes_out);
    file_offset = 0;
  }
  head_count = (size_t) (headc + save_num_procs) + (size_t) num_trees;
  zpadding = (align - (head_count * sizeof (uint64_t)) % align) % align;
  if (zpadding > 0 || rank > 0) {
    retval = sc_io_source_read (src, NULL, (long)
                                (file_offset + zpadding +
                                 gfq[rank] * comb_size), NULL);
    SC_CHECK_ABORT (!retval, "seek data");
  }

  /* read quadrant coordinates and data interleaved */
  qarr =
    sc_array_new_size (sizeof (p4est_qcoord_t), (P4EST_DIM + 1) * zcount);
  qap = (p4est_qcoord_t *) qarr->array;
  darr = NULL;
  dap = NULL;
  lbuf = NULL;
  if (load_data) {
    P4EST_ASSERT (data_size == save_data_size && data_size > 0);
    darr = sc_array_new_size (data_size, zcount);
    dap = darr->array;
    lbuf = P4EST_ALLOC (char, comb_size);
  }
  for (zz = 0; zz < zcount; ++zz) {
    if (load_data) {
      retval = sc_io_source_read (src, lbuf, comb_size, NULL);
      SC_CHECK_ABORT (!retval, "read quadrant with data");
      memcpy (qap, lbuf, qbuf_size);
      memcpy (dap, lbuf + qbuf_size, data_size);
    }
    else {
      retval = sc_io_source_read (src, qap, qbuf_size, NULL);
      SC_CHECK_ABORT (!retval, "read quadrant with data");
      if (save_data_size > 0) {
        retval = sc_io_source_read (src, NULL, save_data_size, NULL);
        SC_CHECK_ABORT (!retval, "seek over data");
      }
    }
    qap += P4EST_DIM + 1;
    dap += data_size;
  }
  P4EST_FREE (lbuf);

  /* seek every process to the end of the source (in case there is data
   * appended to the end of this source) */
  if (gfq[num_procs] > gfq[rank + 1]) {
    retval = sc_io_source_read
      (src, NULL, (long) (gfq[num_procs] - gfq[rank + 1]) * comb_size, NULL);
    SC_CHECK_ABORT (!retval, "seek to end of data");
  }

  /* create p4est from accumulated information */
  p4est = p4est_inflate (mpicomm, conn, gfq, pertree,
                         qarr, darr, user_pointer);
  sc_array_destroy (qarr);
  if (darr != NULL) {
    sc_array_destroy (darr);
  }
  P4EST_FREE (pertree);
  P4EST_FREE (gfq);

  /* assert that we loaded a valid forest and return */
  SC_CHECK_ABORT (p4est_is_valid (p4est), "invalid forest");

  return p4est;
}
