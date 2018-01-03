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

/*
 * This is an internal header file.
 * It is included from a c file that may be preceded by p4est_to_p8est.h.
 */

#ifndef PART_GLOBAL_H
#define PART_GLOBAL_H

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif /* P4_TO_P8 */

SC_EXTERN_C_BEGIN;

typedef double      (*part_init_density_t) (double x, double y, double z,
                                            void *data);

typedef struct part_global
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  int                 minlevel;
  int                 maxlevel;
  int                 bricklev;
  int                 order;
  int                 olap_notify;
  int                 ntop;
  int                 nint;
  int                 nbot;
  int                 vtk;
  int                 checkp;
  int                 printn;
  double              num_particles;
  double              elem_particles;
  double              deltat;
  double              finaltime;
  const char         *prefix;

  int                 bricklength;
  int                 stage;
  double              global_density;
  double              t, lxyz[3], hxyz[3], dxyz[3];
  p4est_locidx_t      prevlp, prev2;
  sc_array_t         *padata;   /**< pa_data_t Numerical data of local particles */
  sc_array_t         *pfound;   /**< int Target rank for local particles */
  sc_array_t         *iremain;  /**< locidx_t Index into padata of stay-local particles */
  p4est_locidx_t      ireindex, ire2;   /**< Running index into iremain */
  p4est_locidx_t      qremain;  /**< Number of particles remaining in quadrant */
  sc_array_t         *ireceive; /**< Index into particle receive buffer */
  p4est_locidx_t      irvindex, irv2;   /**< Running index into ireceive */
  p4est_locidx_t      qreceive; /**< Number of particles received for quadrant */
  sc_array_t         *recevs;   /**< comm_prank_t with one entry per receiver, sorted */
  sc_array_t         *recv_req; /**< sc_MPI_Request for receiving */
  sc_array_t         *send_req; /**< sc_MPI_Request for sending */
  sc_array_t         *prebuf;   /**< pa_data_t All received particles */
  sc_array_t         *cfound;   /**< char Flag for received particles */
  sc_hash_t          *psend;    /**< comm_psend_t with one entry per receiver */
  sc_mempool_t       *psmem;    /**< comm_psend_t to use as hash table entries */
  p4est_locidx_t      qcount;   /**< Count local quadrants in partition callback */
  sc_array_t         *src_fixed;        /**< int Particle counts per quadrant */
  sc_array_t         *dest_fixed;       /**< int Particle counts per quadrant */
  part_init_density_t pidense;
  void               *piddata;
  p4est_gloidx_t      gpnum, gplost;

  p4est_connectivity_t *conn;
  p4est_t            *p4est;
}
part_global_t;

SC_EXTERN_C_END;

#endif /* !PART_GLOBAL_H */
