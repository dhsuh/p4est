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

#if !defined(P4EST_NEIGH_H)
#define      P4EST_NEIGH_H

#include <p4est.h>

/** An object that can be used to manage neighborhood collective
 * communication.  A p4est_neigh_t can be used to do some of the
 * MPI_Neighborhood_* operations (in the future, maybe all of them).  Ideally,
 * MPI implementations of MPI Neighborhood collectives will be flexible in
 * their algorithms that \a p4est_neigh_t becomes redundant.  But as of this
 * writing in 2018, popular MPI implementations only have loop-over-neighbors
 * implementations of neighborhoods */
typedef struct p4est_neigh_s p4est_neigh_t;

/** Create a neighborhood. The communicator from the p4est will be duplicated so
 * that the neighborhood can manage its own tags without conflict.  The reason
 * that a p4est is passed instead of just a communicator is that the
 * information about the subsets of the space-filling curve owned by each
 * process can be used to scheduling the communication pattern of neighborhood
 * collectives.
 *
 * \param[in] p4est A forest.
 * \param[in] n_neigh number of neighbors.
 * \param[in] neigh_procs the list of neighbors.  Will be copied, the user can
 *            modify or delete this list after the call to
 *            \a p4est_neigh_new().
 *
 * \return A p4est neighborhood.
 */
p4est_neigh_t *p4est_neigh_new (p4est_t *p4est, int n_neigh, const int *neigh_procs);

/** Destroy a neighborhood object.
 *
 * \param[in] neigh  The neighborhodd object to be destroyed.
 */
void p4est_neigh_destroy (p4est_neigh_t *neigh);

/** Get the array of procs.
 *
 * \param [in] neigh 	Neighborhood object from which array of neighboring 
 *	                    processors will be fetched from.
 */
void p4est_neigh_get_procs (p4est_neigh_t *neigh, int *n_neigh,
                            const int **neigh_procs);

void p4est_neigh_allgather (p4est_neigh_t *neigh,
                            sc_array_t *send_array,
                            sc_array_t *recv_array,
                            int ordered);

void p4est_neigh_allgatherv (p4est_neigh_t *neigh,
                             sc_array_t *send_array,
                             sc_array_t *recv_buf,
                             sc_array_t *recv_offsets);

void p4est_neigh_alltoall (p4est_neigh_t *neigh,
                           sc_array_t *send_array,
                           sc_array_t *recv_array,
                           int ordered);

void p4est_neigh_alltoallv (p4est_neigh_t *neigh,
                            sc_array_t *send_array,
                            sc_array_t *send_offsets,
                            sc_array_t *recv_buf,
                            sc_array_t *recv_offsets);



#endif /* P4EST_NEIGH_H */

