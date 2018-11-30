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

#if !defined(P4_TO_P8)
#include <p4est_neigh.h>
#else
#include <p8est_neigh.h>
#endif

typedef enum
{
  P4EST_NEIGH_DEFAULT = -1,
  P4EST_NEIGH_BASIC = 0,
  P4EST_NEIGH_MPI,
  P4EST_NEIGH_TREE,
  P4EST_NEIGH_NUM_METHODS,
}
p4est_neigh_method_t;

enum
{
  P4EST_NEIGH_BASIC_ALL,
  P4EST_NEIGH_BASIC_ALLV,
};

struct p4est_neigh_s
{
  sc_MPI_Comm comm;
  int size;
  int rank;
  int n_neigh;
  int *neigh_procs;
  int *neigh_procs_sorted;
  int *neigh_perm;
  sc_array_t *reqs;
  sc_array_t *msgs;
  p4est_neigh_method_t method;
};

typedef struct p4est_neigh_msg_s
{
  int         source;
  size_t      count;
#if defined(P4EST_ENABLE_MPIMPROBE)
  MPI_Message msg;
#endif
}
p4est_neigh_msg_t;

static int
p4est_neigh_msg_cmp (const void *A, const void *B)
{
  const p4est_neigh_msg_t *msgA = (const p4est_neigh_msg_t *) A;
  const p4est_neigh_msg_t *msgB = (const p4est_neigh_msg_t *) B;

  return (msgA->source - msgB->source);
}

p4est_neigh_t *
p4est_neigh_new (p4est_t *p4est, int n_neigh, const int *neigh_procs)
{
  int mpiret;

  p4est_neigh_t *neigh;

  neigh = P4EST_ALLOC_ZERO (p4est_neigh_t, 1);
  neigh->method = P4EST_NEIGH_BASIC;
  
  P4EST_ASSERT (n_neigh >= 0);

#if !defined(P4EST_ENABLE_MPI)
  P4EST_ASSERT (n_neigh == 0 || n_neigh == 1);
  P4EST_ASSERT (n_neigh == 0 || neigh_procs[0] == 0);
#endif

  if (neigh->method == P4EST_NEIGH_MPI) {
    /* TODO: create distributed graph MPI_Dist_graph_create_adjacent */
	mpiret = MPI_Dist_graph_create_adjacent(p4est->mpicomm, n_neigh,
											neigh_procs, MPI_UNWEIGHTED,
											n_neigh, neigh_procs, 
											MPI_UNWEIGHTED, NULL, 0, &(neigh->comm));
	SC_CHECK_MPI (mpiret);
  }
  else {
    mpiret = sc_MPI_Comm_dup (p4est->mpicomm, &(neigh->comm));
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Comm_size (neigh->comm, &(neigh->size));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (neigh->comm, &(neigh->rank));
  SC_CHECK_MPI (mpiret);

  neigh->n_neigh = n_neigh;
  neigh->neigh_procs = P4EST_ALLOC (int, n_neigh);
  memcpy (neigh->neigh_procs, neigh_procs, n_neigh);

  {
    sc_array_t neigh_array;
    int        is_sorted;

    sc_array_init_data (&neigh_array, neigh->neigh_procs, sizeof(int),
                        n_neigh);

    is_sorted = sc_array_is_sorted (&neigh_array, sc_int_compare);

    if (!is_sorted) {
      int *neigh_sorter = P4EST_ALLOC (int, 2 * n_neigh), nq;
      neigh->neigh_procs_sorted = P4EST_ALLOC (int, n_neigh);
      neigh->neigh_perm = P4EST_ALLOC (int, n_neigh);

      for (nq = 0; nq < n_neigh; nq++) {
        neigh_sorter[2 * nq] = neigh->neigh_procs[nq];
        neigh_sorter[2 * nq + 1] = nq;
      }
      memcpy (neigh->neigh_procs_sorted, neigh->neigh_procs, n_neigh * sizeof(int));
      sc_array_init_data (&neigh_array, neigh_sorter, 2 * sizeof(int),
                          n_neigh);
      sc_array_sort (&neigh_array, sc_int_compare);
      for (nq = 0; nq < n_neigh; nq++) {
        int orig = neigh_sorter[2 * nq + 1];

        neigh->neigh_procs_sorted[nq] = neigh_sorter[2 * nq];
        neigh->neigh_perm[orig] = nq;
      }
      P4EST_FREE (neigh_sorter);
    }
    else {
      neigh->neigh_procs_sorted = neigh->neigh_procs;
    }
  }

#if defined(P4EST_ENABLE_DEBUG)
  {
    int *neigh_in, *neigh_out, nq, diff;

    neigh_in = P4EST_ALLOC_ZERO(int, neigh->size);
    neigh_out = P4EST_ALLOC(int, neigh->size);

    for (nq = 0; nq < neigh->n_neigh; nq++) {
      neigh_in[nq] = 1;
    }

    mpiret = sc_MPI_Alltoall (neigh_in, 1, MPI_INT, neigh_out, 1, MPI_INT,
                              neigh->comm);
    SC_CHECK_MPI (mpiret);

    diff = memcmp (neigh_in, neigh_out, neigh->n_neigh * sizeof(int));

    P4EST_ASSERT (!diff);
  }
#endif

  neigh->reqs = sc_array_new (sizeof (sc_MPI_Request));
  neigh->msgs = sc_array_new (sizeof (p4est_neigh_msg_t));

  return neigh;
}

void
p4est_neigh_destroy (p4est_neigh_t *neigh)
{
  int mpiret;

  sc_array_destroy (neigh->msgs);
  sc_array_destroy (neigh->reqs);
  mpiret = sc_MPI_Comm_free(&(neigh->comm));
  SC_CHECK_MPI (mpiret);
  P4EST_FREE (neigh->neigh_perm);
  if (neigh->neigh_procs_sorted != neigh->neigh_procs) {
    P4EST_FREE (neigh->neigh_procs_sorted);
  }
  P4EST_FREE (neigh->neigh_procs);

  P4EST_FREE (neigh);
}

void
p4est_neigh_get_procs (p4est_neigh_t *neigh,
                       int *n_neigh,
                       const int **neigh_procs)
{
  if (n_neigh) {
    *n_neigh = neigh->n_neigh;
  }
  if (neigh_procs) {
    *neigh_procs = neigh->neigh_procs;
  }
}

static void
p4est_neigh_all_basic (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *recv_array,
                       int ordered, int alltoall)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs;
  int nq;
  int mpiret;

  reqs = (sc_MPI_Request *) sc_array_push_count (neigh->reqs, nreqs);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  for (nq = 0; nq < n_neigh; nq++) {
    size_t idx = alltoall ? nq: 0;
    mpiret = sc_MPI_Isend (sc_array_index (send_array, idx),
                           (int) send_array->elem_size,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALL,
                           comm,
                           &reqs[2 * nq]);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Irecv (sc_array_index (recv_array, nq),
                           (int) recv_array->elem_size,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALL,
                           comm,
                           &reqs[2 * nq + 1]);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Waitall (nreqs, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* release requests */
  neigh->reqs->elem_count -= nreqs;
}

void
p4est_neigh_allgather (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *recv_array,
                       int ordered)
{
  sc_array_t send_array_view, recv_array_view;
  size_t total_size;

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == neigh->n_neigh * send_array->elem_count);

  total_size = send_array->elem_size * send_array->elem_count;
  if (!total_size) {
    return;
  }

  if (send_array->elem_count > 1) {
    sc_array_init_data (&send_array_view, send_array->array, send_array->elem_size * send_array->elem_count, 1);
    sc_array_init_data (&recv_array_view, recv_array->array, send_array->elem_size * send_array->elem_count, neigh->n_neigh);
    send_array = &send_array_view;
    recv_array = &recv_array_view;
  }

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (send_array->elem_count == 1);
  P4EST_ASSERT (recv_array->elem_count == neigh->n_neigh);

  p4est_neigh_all_basic (neigh, send_array, recv_array, ordered, 0);
}

#if defined(P4EST_ENABLE_MPIMPROBE)
static void
p4est_neigh_allv_basic_mprobe (p4est_neigh_t *neigh,
                               sc_array_t *send_array,
                               sc_array_t *send_offsets,
                               sc_array_t *recv_array,
                               sc_array_t *recv_offsets,
                               int alltoall)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nmsgs = n_neigh;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs;
  size_t *offsets = NULL;
  p4est_neigh_msg_t *msgs;
  int nq;
  int mpiret;
  int sorted = (recv_offsets != NULL);
  size_t off;
  size_t send_count = 0, recv_count;
  size_t elem_size = send_array->elem_size;
  char *send_buf = NULL;
  char *recv_buf = NULL;
  int *neigh_perm = neigh->neigh_perm;

  reqs = (sc_MPI_Request *) sc_array_push_count (neigh->reqs, nreqs);
  msgs = (p4est_neigh_msg_t *) sc_array_push_count (neigh->msgs, nmsgs);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  if (alltoall) {
    offsets = sc_array_index (send_offsets, 0);
    send_buf = send_array->elem_count ? sc_array_index (send_array, 0) : NULL;
  }
  if (n_neigh && send_array->elem_count && !alltoall) {
    send_count = send_array->elem_count;
    send_buf = sc_array_index (send_array, 0);
  }
  /* send messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *buf = alltoall ? &send_buf[offsets[nq] * elem_size] : send_buf;
    size_t count = alltoall ? (offsets[nq + 1] - offsets[nq]) * elem_size : send_count * elem_size;

    mpiret = sc_MPI_Isend (buf,
                           count,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALLV,
                           comm,
                           &reqs[2 * nq]);
    SC_CHECK_MPI (mpiret);
  }

  offsets = recv_offsets ? (size_t *) sc_array_index (recv_offsets, 0) : NULL;
  /* get sizes with mprobe */
  for (nq = 0, recv_count = 0; nq < n_neigh; nq++) {
    MPI_Status status;
    int        count;

    mpiret = MPI_Mprobe (sc_MPI_ANY_SOURCE,
                         P4EST_NEIGH_BASIC_ALLV,
                         comm,
                         &(msgs[nq].msg),
                         &status);
    SC_CHECK_MPI (mpiret);
    msgs[nq].source = status.MPI_SOURCE;
    mpiret = MPI_Get_count (&status, MPI_BYTE, &count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (count % elem_size == 0);
    count = count / elem_size;
    msgs[nq].count = (size_t) count;
    recv_count += (size_t) count;
  }
  off = recv_array->elem_count;
  (void) sc_array_push_count (recv_array, recv_count);
  if (recv_array->elem_count) {
    recv_buf = recv_array->array;
  }
  if (sorted) {
    sc_array_t msg_sorter;

    sc_array_init_data (&msg_sorter, msgs, sizeof (p4est_neigh_msg_t),
                        nmsgs);
    sc_array_sort (neigh->msgs, p4est_neigh_msg_cmp);
  }
  for (nq = 0; nq < n_neigh; nq++) {
    int q = (sorted && neigh_perm) ? neigh_perm[nq] : nq;
    p4est_neigh_msg_t *msg = &msgs[q];

    P4EST_ASSERT (msg->source == neigh->neigh_procs[nq]);
    mpiret = MPI_Imrecv (&recv_buf[off * elem_size],
                         msg->count * elem_size, sc_MPI_BYTE, &(msg->msg),
                         &reqs[2 * nq + 1]);
    SC_CHECK_MPI (mpiret);
    if (offsets) {
      offsets[nq] = off;
    }
    off += msg->count;
  }
  if (offsets) {
    offsets[n_neigh] = off;
  }
  P4EST_ASSERT (off == recv_array->elem_count);
  mpiret = sc_MPI_Waitall (nreqs, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* release requests */
  if (!sorted) {
    recv_offsets->elem_count -= n_neigh + 1;
  }
  neigh->msgs->elem_count -= nmsgs;
  neigh->reqs->elem_count -= nreqs;
}
#else
static void
p4est_neigh_allv_basic_counts (p4est_neigh_t *neigh,
                               sc_array_t *send_array,
                               sc_array_t *send_offsets,
                               sc_array_t *recv_array,
                               sc_array_t *recv_offsets,
                               int alltoall)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs;
  size_t *soffsets = NULL;
  size_t *roffsets = NULL;
  int nq;
  int mpiret;
  int sorted = (recv_offsets != NULL);
  size_t off;
  size_t send_count = 0, recv_count;
  size_t elem_size = send_array->elem_size;
  char *send_buf = NULL;
  char *recv_buf = NULL;
  sc_array_t *send_counts, *recv_counts;

  send_counts = sc_array_new_size (sizeof (size_t), alltoall ? neigh->n_neigh : 1);
  recv_counts = sc_array_new (sizeof (size_t));
  if (alltoall) {
    for (nq = 0; nq < n_neigh; nq++) {
      size_t offstart = *((size_t *) sc_array_index (send_offsets, nq));
      size_t offend = *((size_t *) sc_array_index (send_offsets, nq+1));

      *((size_t *) sc_array_index (send_counts, nq)) = offend - offstart;
    }
  }
  else {
    *((size_t *) sc_array_index (send_counts, 0)) = send_array->elem_count;
  }

  p4est_neigh_all_basic (neigh, send_counts, recv_counts, 1, alltoall);

  for (nq = 0, recv_count = 0; nq < n_neigh; nq++) {
    recv_count += *((size_t *) sc_array_index (recv_counts, nq));
  }

  reqs = (sc_MPI_Request *) sc_array_push_count (neigh->reqs, nreqs);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  if (alltoall) {
    soffsets = sc_array_index (send_offsets, 0);
  }
  if (sorted) {
    roffsets = sc_array_index (recv_offsets, 0);
  }
  if (send_array->elem_count) {
    send_buf = send_array->array;
  }
  off = recv_array->elem_count;
  (void) sc_array_push_count (recv_array, recv_count);
  if (recv_array->elem_count) {
    recv_buf = recv_array->array;
  }
  /* send / recv messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *sbuf = alltoall ? &send_buf[soffsets[nq] * elem_size] : send_buf;
    size_t scount = alltoall ? (soffsets[nq + 1] - soffsets[nq]) * elem_size : send_count * elem_size;
    char *rbuf = &recv_buf[off * elem_size];
    size_t rcount = *((size_t *) sc_array_index (recv_counts, nq));

    if (scount) {
      mpiret = sc_MPI_Isend (sbuf,
                             scount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq]);
      SC_CHECK_MPI (mpiret);
    }
    if (rcount) {
      mpiret = sc_MPI_Irecv (rbuf,
                             rcount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq + 1]);
      SC_CHECK_MPI (mpiret);

      if (recv_offsets) {
        roffsets[nq] = off;
      }
    }
    off += rcount;
  }
  if (recv_offsets) {
    roffsets[nq] = off;
  }
  P4EST_ASSERT (off == recv_array->elem_count);

  sc_array_destroy (recv_counts);
  sc_array_destroy (send_counts);
}
#endif

static void
p4est_neigh_allv_basic (p4est_neigh_t *neigh,
                        sc_array_t *send_array,
                        sc_array_t *send_offsets,
                        sc_array_t *recv_array,
                        sc_array_t *recv_offsets,
                        int alltoall)
{
#if defined(P4EST_ENABLE_MPIMPROBE)
  p4est_neigh_allv_basic_mprobe (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall);
#else
  p4est_neigh_allv_basic_counts (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall);
#endif
}

void
p4est_neigh_allgatherv (p4est_neigh_t *neigh,
                        sc_array_t *send_array,
                        sc_array_t *recv_array,
                        sc_array_t *recv_offsets)
{
  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (recv_array));
  P4EST_ASSERT (recv_offsets == NULL || recv_offsets->elem_size == sizeof(size_t));
  P4EST_ASSERT (recv_offsets == NULL || recv_offsets->elem_count == neigh->n_neigh + 1);

  if (!send_array->elem_size) {
    return;
  }

  p4est_neigh_allv_basic (neigh, send_array, NULL, recv_array, recv_offsets, 0);
}

void
p4est_neigh_alltoall (p4est_neigh_t *neigh,
                      sc_array_t *send_array,
                      sc_array_t *recv_array,
                      int ordered)
{
  sc_array_t send_array_view, recv_array_view;

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == send_array->elem_count);
  P4EST_ASSERT (send_array->elem_count % neigh->n_neigh == 0);

  if (!send_array->elem_size) {
    return;
  }

  if (send_array->elem_count > neigh->n_neigh) {
    size_t mult = send_array->elem_count / neigh->n_neigh;

    sc_array_init_data (&send_array_view, send_array->array, send_array->elem_size * mult, neigh->n_neigh);
    sc_array_init_data (&recv_array_view, recv_array->array, send_array->elem_size * mult, neigh->n_neigh);
    send_array = &send_array_view;
    recv_array = &recv_array_view;
  }

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == send_array->elem_count);
  P4EST_ASSERT (send_array->elem_count == neigh->n_neigh);

  p4est_neigh_all_basic (neigh, send_array, recv_array, ordered, 1);
}

void
p4est_neigh_alltoallv (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *send_offsets,
                       sc_array_t *recv_array,
                       sc_array_t *recv_offsets)
{
  int mpiret;
  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (recv_array));
  P4EST_ASSERT (send_offsets->elem_size == sizeof (size_t));
  P4EST_ASSERT (send_offsets->elem_count == neigh->n_neigh + 1);
  P4EST_ASSERT (recv_offsets == NULL || recv_offsets->elem_count == neigh->n_neigh + 1);

  if (!send_array->elem_size) {
    return;
  }

  if (neigh->method == P4EST_NEIGH_MPI) {
      mpiret = MPI_Neighbor_alltoallv (send_array->array,   send_array->elem_count, 
	    	  						   send_offsets->array, sc_MPI_BYTE, 
			    					   recv_array->array,   recv_array->elem_count,
				    				   recv_offsets->array, sc_MPI_BYTE, neigh->comm);
      }
  else {
      p4est_neigh_allv_basic (neigh, send_array, send_offsets, recv_array, recv_offsets, 1);
  }
}
