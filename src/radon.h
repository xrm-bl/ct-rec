/* radon.h : forward Radon transform of one CT slice (rec2rec).
 *
 * img  : N x N float slice (row-major, v*N+h)
 * sino : M x N float sinogram (row m = view, column r = detector bin)
 * M    : number of views over 180 deg; view m is at angle t0 + pi*m/M
 * t0   : offset angle [rad] (RA0)
 *
 * Geometry matches the CBP backprojector (cbp_thread.c / cbp.cu):
 *   r(h,v) = (h-N12)*cos(th) - (v-N12)*sin(th) + N12,   N12=(N-1)/2
 * i.e. the rotation axis is the image centre and the sinogram is to be
 * reconstructed with CBP(1.0, -(N-1)/2, t0).
 *
 * Implementations: radon_omp.c (CPU, OpenMP) / radon_g.cu (GPU kernel).
 */
extern void	RadonSlice(const float *img,int N,float *sino,int M,double t0);
