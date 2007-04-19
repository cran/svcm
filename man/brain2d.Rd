\name{brain2d}
\alias{brain2d}
\docType{data}
\title{Two-dimensional Diffusion Weighted Dataset}
\description{
  The data set consists of six transformed diffusion weighted images
  (DWI) showing a representative axial slice of the human brain. The
  stored values can directly be passed to estimate the diffusion tensor
  elements (regression coefficients) using a transform of the applied
  gradients as regressors. 
}
\usage{data(brain2d)}
\format{
  The first two dimensions provide the transformed signal intensities of
  one brain slice sized 90 \eqn{\times}{x} 75 voxels. The third
  dimension encodes for the direction of the six applied diffusion
  weighting gradients.
}
\details{
  The present DTI data set was acquired at 1.5 T (Signa Echospeed; GE
  Medical Systems) using a spin-echo echo-planar sequence with TR/TE =
  4200ms/120ms and diffusion gradients in a six noncollinear directions
  with a b-value of 880 s/mm\eqn{\mbox{\textasciicircum}2}{^2}. One axial
  slice was selected from a volume of six DWI (b = 880
  s/mm\eqn{\mbox{\textasciicircum}}{^}2) and one reference image (b = 0
  s/mm\eqn{\mbox{\textasciicircum}}{^}2). In-plane resolution amounts to
  1.875 \eqn{\times}{x} 1.875 mm\eqn{\mbox{\textasciicircum}}{^}2.  
  
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
data(brain2d)
dim(brain2d)
old.par <- par(no.readonly=TRUE)
par(pin=c(6, 1.2))
image(t(matrix(brain2d, dim(brain2d)[1], dim(brain2d)[2]*6)), axes=FALSE,
      col=grey.colors(256))
title("Six Diffusion Weighted Images")
par(old.par)  
}
\keyword{datasets}
