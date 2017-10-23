 void unpack_planes(lb_t * lb, int* nlocal, double* fptr );
 void unpack_edges(lb_t * lb, int* nlocal, double* fptr );
 void unpack_corners(lb_t * lb, int* nlocal, double* fptr );


void unpack_planes_i(lb_t * lb, int* nlocal, double* fptr,int i );
void unpack_edges_i(lb_t * lb, int* nlocal, double* fptr, int i );
void unpack_corners_i(lb_t * lb, int* nlocal, double* fptr, int i );

void free_planes(lb_t * lb);
void free_planes_i(lb_t * lb, int i);
void free_edges(lb_t * lb);
void free_edges_i(lb_t * lb, int i);
void free_corners(lb_t * lb);
void free_corners_i(lb_t * lb, int i);
