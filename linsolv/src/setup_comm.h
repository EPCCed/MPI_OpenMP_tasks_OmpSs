/*
 * setup_comm.h
 */

#ifndef SETUP_COMM_H_
#define SETUP_COMM_H_

typedef struct
{
  int thisdomain; /* global ID of this domain */
  int ndomains; /* overall number of domains */
  int ncommdomains; /* number of communication partners */
  int nownpoints; /* number of own points */
  int naddpoints; /* number of additional (halo) points */

  int *addpoint_owner; /* [naddpoints] domain owning halo point */
  int *addpoint_id; /* [naddpoints] local point index in owner domain */
  int *commpartner; /* [ncommdomains] list of communication partners */
  int *sendcount; /* [ndomains] number of pnts to send to domain */
  int *recvcount; /* [ndomains] number of pnts to receive from domain */
  int **sendindex; /* [ndomains][sendcount[k]] indices to send to domain */
  int **recvindex; /* [ndomains][recvcount[k]] indices to recv from domain */

} CommMap ;

/*******************************************************************************
*
*******************************************************************************/
CommMap *init_comap(int thisdom, int ndom);

/*******************************************************************************
*
*******************************************************************************/
void free_comap(CommMap **comap);

/*******************************************************************************
*
*******************************************************************************/
void create_recvsend_index(CommMap *comap);

/*******************************************************************************
*
*******************************************************************************/
void mpi_parallel_init(int argc, char *argv[], int *rank, int *ndom);

/*******************************************************************************
*
*******************************************************************************/
void mpi_parallel_end(void);

#endif /* SETUP_COMM_H_ */
