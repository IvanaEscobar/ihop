!BOP
!    !ROUTINE: SIZE.h
!    !INTERFACE:
!    include SIZE.h
!    !DESCRIPTION: \bv
!     *==========================================================*
!     | SIZE.h Declare size of underlying computational grid.
!     *==========================================================*
!     | The design here supports a three-dimensional model grid
!     | with indices I,J and K. The three-dimensional domain
!     | is comprised of nPx*nSx blocks (or tiles) of size sNx
!     | along the first (left-most index) axis, nPy*nSy blocks
!     | of size sNy along the second axis and one block of size
!     | Nr along the vertical (third) axis.
!     | Blocks/tiles have overlap regions of size OLx and OLy
!     | along the dimensions that are subdivided.
!     *==========================================================*
!     \ev
!
!     Voodoo numbers controlling data layout:
!     sNx :: Number of X points in tile.
!     sNy :: Number of Y points in tile.
!     OLx :: Tile overlap extent in X.
!     OLy :: Tile overlap extent in Y.
!     nSx :: Number of tiles per process in X.
!     nSy :: Number of tiles per process in Y.
!     nPx :: Number of processes to use in X.
!     nPy :: Number of processes to use in Y.
!     Nx  :: Number of points in X for the full domain.
!     Ny  :: Number of points in Y for the full domain.
!     Nr  :: Number of points in vertical direction.
!EOP
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr
      PARAMETER (                                                                                                                       &
     &           sNx =  62,                                                                                                             &
     &           sNy =  62,                                                                                                             &
     &           OLx =   2,                                                                                                             &
     &           OLy =   2,                                                                                                             &
     &           nSx =   1,                                                                                                             &
     &           nSy =   1,                                                                                                             &
     &           nPx =   1,                                                                                                             &
     &           nPy =   1,                                                                                                             &
     &           Nx  = sNx*nSx*nPx,                                                                                                     &
     &           Ny  = sNy*nSy*nPy,                                                                                                     &
     &           Nr  =   15)

!     MAX_OLX :: Set to the maximum overlap region size of any array
!     MAX_OLY    that will be exchanged. Controls the sizing of exch
!                routine buffers.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,                                                                                                    &
     &            MAX_OLY = OLy )
