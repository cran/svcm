\name{brain3d}
\alias{brain3d}
\docType{data}
\title{Three-dimensional Diffusion Weighted Dataset}
\description{
  To keep the computational effort low a volume of 15 \eqn{\times}{x} 30
  \eqn{\times}{x} 6 voxels was chosen from the original whole brain
  volume. The extract depicts the posterior part of the lateral
  ventricles and the corpus callosum so both areas with low and high
  signal intensities are contained. The six transformed diffusion weighted
  images can directly be passed to estimate the diffusion tensor
  elements (regression coefficients) using a transform of the applied
  gradients as regressors.
}
\usage{data(brain3d)}
\format{
  The first three dimensions of the data array contain the number of voxels in
  x-, y- and z-direction. The fourth dimension encodes for the direction
  of the six applied diffusion weighting gradients. 
}
\details{
  The present DTI data set was acquired at 1.5 T (Signa Echospeed; GE
  Medical Systems) using a spin-echo echo-planar sequence with TR/TE =
  4200ms/120ms and diffusion gradients in a six noncollinear directions
  with a b-value of 880 s/mm\eqn{\mbox{\textasciicircum}}{^}2. The
  extracted volume originates from six DWI (b = 880
  s/mm\eqn{\mbox{\textasciicircum}}{^}2) and one reference image (b = 0
  s/mm\eqn{\mbox{\textasciicircum}}{^}2). In-plane resolution amounts to
  1.875 \eqn{\times}{x} 1.875 mm\eqn{\mbox{\textasciicircum}}{^}2, slice
  thickness is 4.0 mm.

  The transformation of the raw signal intensities, \deqn{ y =
    -\frac{1}{b}\log\left(\frac{S_i}{S_0}\right), i=1,\ldots,6}{y = - 1/b
    log(S_i/S_0)} is derived from the Stejskal-Tanner equation and is
  proposed by Papadakis et al.
}
\source{
  Diffusion Tensor Imaging was performed at the Max-Planck-Institute of
  Psychiatry, Munich, Germany.
}
\references{
  Basser P. J. and Jones D. K. (2002) Diffusion-tensor MRI: Theory,
  experimental design and data analysis -- a technical review. \emph{NMR
    in Biomedicine} \bold{15}, 456-467.

  Mori S. and Barker P. B. (1999) Diffusion magnetic resonance imaging:
  Its principle and applications. \emph{The Anatomical Record}
  \bold{257}, 102-109.

  Papadakis N. G., Xing D., Huang C. L.-H., Hall L. and Carpenter
  T. A. (1999). A comparative study of acquisition schemes for diffusion
  tensor imaging using MRI. Journal of Magnetic Resonance 137, 67-82. 
  
  Stejskal E. O. and Tanner J. E. (1965) Spin diffusion measurements:
  Spin echoes in the presence of time-dependent field
  gradient. \emph{The Journal of Chemical Physics} \bold{42}, 288-292.
  
}
\examples{
data(brain3d)
dim(brain3d)
old.par <- par(no.readonly = TRUE) 
par(pin=c(1.1, 3.4), mfrow=c(1, 6))
for (i in 1:dim(brain3d)[4])
  image(matrix(aperm(brain3d[,,,i], c(2,1,3)), nrow=dim(brain3d)[2]),
        axes=FALSE, col=grey.colors(256), main=paste("DWI", i))
title("6 DWIs of a (15 x 30 x 6) human brain extract in axial view",
      outer=TRUE, line=-10) 
par(old.par)
}
\keyword{datasets}
