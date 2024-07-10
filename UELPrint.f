!Abaqus UEL subroutine for linear single element plain stress
!this code might be used for multiple element  
! For multiple element, change the nelem=1 in UEL subroutine part to the UEL element number from your input file
!do the same change in UMAT subroutine part, nelem=1 to te UMAT number of element from your input file
!This code is modified from a previous code (https://github.com/KianAghani/UEL)
!this code will write two text file which might be helpful for data analysis
!SingleElementUEL.inp file and its abaqus_v6.env file are provied in the directory
!fortran free formet env file 
!it will print several value in console which will help to understand how abaqus caught different subroutine
!this printing will help to compare the differnet results with matlab results
!the matlab code for the same problem is given below
!don't get frightend, half of the code is comments for understanding 


!--------------------Ranning instruction of this file-------------------------------------
!two changes are necessary
!1. filepath parameter (filepath='/home/mmollah/Documents/2024June12/'))in the module WorkDir need to be changed to working directory
!2. filename(filename = '/home/mmollah/Documents/2024June10/UELtoUMATResultsTransfer.txt') in the Subroutine UMAT need to be changed to working directroy
!need to put everything (input file, this fortran code, environment file in the working directory
!ensure that you are at the working directory at your linux or windows terminal
!syntax for ranning the code is: (no file extension is needed in the syntax)
!abaqus job=SingleElementUEL user=UELPrint -inter

!-----------------------------------------------------------------------------------------------------------------
!--- Module for transfering the UEL stress to the dummy element umat
module kvisual
    implicit none
    real*8 sigout(1000,4,4)
    save
end module   
!----------------------------------------------------------------------
!--- Module specifying the output file namne
module WorkDir
    !
    !----Define the file path and root name

    character*80 filepath
    parameter (filepath='/home/mmollah/Documents/2024June12/')

    integer unit_number
    parameter (unit_number=200)

    character*80 OutputFileName
    parameter (OutputFileName='ResultsUEL.dat')
end module
!----------------------------------------------------------------------
    !-- Check aqaqus documentation on uexternaldb https://help.3ds.com/2024x/english/dsdoc/sima3dxsubrefmap/simasub-c-uexternaldb.htm?contextscope=onpremise&highlight=UEXTERNALDB&id=bd9af4e4ac2e48419ffb8a57da307ede&analyticsContext=search-result&analyticsSearch=UEXTERNALDB&myapp=false
    
    !LOP=0 indicates that the user subroutine is being called at the start of the analysis.
    !LOP=4 indicates that the user subroutine is being called at the beginning of a restart analysis. When LOP=4, all necessary external files should be opened and properly positioned and all information required for the restart       
     !should be read from the external files.
!----------------------------------------------------------------------
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!external database management
    use WorkDir
    implicit none

    integer lop, lrestart, kstep, kinc

    real*8 TIME(2), dtime

    !variable type declaration for text file   
    character*80 filename
    integer ierr, length1, length2
    print *, '@@@@@@@@@@@@@@caught of UEXTERNALDB@@@@@@@@@@@@@@'



    ! Only open the file at the start of the analysis
    if (lop==0 .or. lop==4) then 
        ! Combine filePath and OutputFileName to form the full filename
        !call getcwd(filepath)
        length1 = index(filePath,' ') - 1
        length2 = index(OutputFileName,' ') - 1
        filename = filePath(1:length1)//OutputFileName(1:length2)
        write(*,*) 'filename=', filename
	
	! Open the file for writing
	!The variable unit_number should be an integer that uniquely identifies the file within the program. It is used in subsequent I/O operations to refer to this file.
	!The variable filename should be a character string that contains the name of the file. This can be a relative or absolute path.
	!unknown' means that the status of the file is not known. If the file exists, it will be opened. If it does not exist, a new file will be created.
	!sequential' means that the file will be accessed sequentially. This is typical for most text file operations.
        open(unit=unit_number, file=filename, status='unknown',  &
                                     access='sequential')
        !rewind statement in Fortran is used to position(file pointer back to the beginning of the file) a file back to its beginning
        !sequential access.It ensures that subsequent read or write operations start from the beginning of the file.
        rewind(unit_number)
        ! Write header information to the file
        !KINC, TIME(1), JELEM, KINTK, SIG(1), SIG(2), SIG(3), SIG(4)
        write(unit_number,'(A)') 'this results is from directly UEL sig variable'
        write(unit_number,*) 'KINC-----TIME(1)-JELEM-KINTK----Sig11------Sig22-------Sig33-------Sig12'
    endif

    RETURN
END 
!----------------------------------------------------------------------

SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,      &
    PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,   &
    KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF, &
    LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

    use kvisual
    use WorkDir
    INCLUDE 'ABA_PARAM.INC'
    !implicit none

	  
    real*8 RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),        &
    SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),       &
    DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),         &
    JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),         &
    PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


    integer maxDOF, ndof, ndim,ninpt, nelem, nsdv,nsvint,ntens, kintk, kinc, kstep
    PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=4,nelem=1,nsdv=4,nsvint=4,ntens=4)   
  
    real*8 dshape(ndim,4),xjac(maxDOF,maxDOF),xjaci(maxDOF,maxDOF),                  &
    bmat(3,2*nnode),bmatT(2*nnode,3),kmat(2*nnode,2*nnode),C(3,3),                      &
    deps(3),stress(3),SIG(4),SDV(4),statevLocal(nsvint)                                 

    real*8 db(3,2*nnode),XYDerivation(ndim,4),force(8)
	
	!variable type declaration for text file   
	character*80 filename
	integer ierr, k1, k, i, j, jelem

    real*8 E, xnu, xa, xb, xc, djac, shapef(4)
	
	!reading in material parameter and initialize variables
	print *, '*******************start of UEL*****************************'
	rhs(:, 1)=0.d0
	amatrx=0.d0
	C=0.d0					
	force=0.d0

	!Elasticity matrix
	E = props(1)
	xnu = props(2)
	xa = E/(1.d0- xnu**2)
	xb = xnu
	xc = (1.d0-xnu)/2.d0

	C(1, 1) = xa
	C(1, 2) = xa*xb
	C(2, 1) = C(1, 2)
	C(2, 2) = C(1, 1)
	C(3, 3) = xa*xc
	print *, 'Material Stifness matirx'
	write(*,*) C
	 
    !looping over total number of integration points, ninpt=4, this is each element loop
    do kintk = 1, ninpt	 
    		print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%start of UEL do loop%%%%%%%%%%%%%%%%%%%%%%'
                deps=0.	           !used in strain increment
		statevLocal=0.d0
		call shapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
		
    	        djac = 1.d0 !determinant of jacobian is initialized as 1 NOT 0(WHY? MATH ERROR) and double precision
		call jacobian(ndim,COORDS,nnode,djac,dshape,xjaci)
	    
		call Bmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)
		
		!SVARS is a global array that holds the state variables for all integration points of all elements.
		!SVARS maintains the state of the entire model throughout the analysis, ensuring continuity between increments and steps.
		!statevLocal is a local array that holds the state variables for the current integration point of the current element.
		!It temporarily stores state variables during computations within the UEL subroutine for a specific integration point.
		
		! Copy state variables from global to local array for the corresponding element, icopy = 1
		!copies the state variables for the current integration point from the global array (SVARS) to the local array (statevLocal).
		
	    	call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,1)
		
		!initialization for variable
		SIG=0.d0
		stress=0.d0
		
		!global to local tranfer for state variable for the corresponding element
		do k1=1,nsdv     !nsdv=4
			SIG(k1) = statevLocal(k1) 
		end do
		
		print *, 'Stress Initialization'
		print *, 'Stress 1', stress(1)
		print *, 'Stress 2', stress(2)
		print *, 'Stress 3', stress(3)
		
		!The strain increment deps is computed using the straininc subroutine. 
		!The stress vector is updated based on the strain increment.  
		call straininc(ntens,ndof,ndim,nnode,MLVARX,bmat,DU,deps)		
	 	
	 	!C=(Stress/Strain) here Strain=deps
	 	!this stress will be used in SIG
		stress=stress+matmul(C,deps)
		
		!initialization
		force=0.d0
		!need to see the documentations
		force=matmul(bmatT,stress)
		
		!The internal force vector RHS is computed by multiplying the transposed strain-displacement matrix with the stress vector.
		do k1=1,8
			!do k2=1, 3
				RHS(k1, 1) = RHS(k1, 1) - force(k1)
			!end do
		end do
		
		!this db will be used in AMATRX. It will simplify the AMATRX calculation
		db=0.
	   	do i=1,3
			do j=1,8
				do k=1,3
					db(i,j)=db(i,j)+C(i,k)*bmat(k,j)*djac
				end do	
			end do
		end do
		
	   	!Compute Element Stiffness Matrix (AMATRX)
		do i=1,8
			do j=1,8
				do k=1,3
					AMATRX(i,j)=AMATRX(i,j)+bmatT(i,k)*db(k,j)
				end do	
			end do
		end do
		!lets print the AMATRX to the console
		!print *, 'Abaqus default results for AMATRX!'
		!write(*,*) AMATRX
		print *, 'Original 8*8 formet of AMATRX'
		write (*,"(8(1x,F16.3))") transpose(AMATRX)	 !this syntax will print in in origina 8*8 formet
		
		!stress is converting into another variable and initialization of that variable
		SIG(1)=stress(1)  !S11
		SIG(2)=stress(2)  !S22
		SIG(3)=0.d0       !S33 in 2D no S33 value
		SIG(4)=stress(3)  !S12
		
		!putting stress value(SIG) of each element to the local state array=statevLocal
	    	do k1=1,nsdv
			statevLocal(k1) = SIG(k1)
		end do
		
		! Copy updated local state variables back to global, icopy = 0
		! copies the updated state variables from the local array (statevLocal) back to the global array (SVARS).
		call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,0)
                
                !putting that global stress to global variable sigout in element and integration point sequential distribution
		do k1=1,nsdv
			sigout(jelem,kintk,k1)=SIG(k1)
		end do
		
		print *, '--SIG from UEL--'
		write (*,*) SIG
		
                !this will write the text file
		write(unit_number, '(I3, F16.8, I3, I3,  4(1x, e12.5))') KINC, TIME(1), JELEM, KINTK, SIG(1), SIG(2), SIG(3), SIG(4)
		print *, '%%%%%%%%%%%%%%%%%%this is the end of UEL do loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
		
	end do       
	print *, '*********************************end of UEL*****************'
	RETURN
END
	  
!******************************************************
!** shape functions                                   *
!******************************************************	
		subroutine shapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
		
		include 'aba_param.inc'
                !implicit none
		
		parameter (gaussCoord=0.577350269189626) 
		real*8 shapef(4),dshape(ndim,4),coord24(2,4)  !ndim=2
		
		!The shapeFunc subroutine computes the shape functions and their derivatives with respect to the natural coordinates for a given integration point
		!kintk: The index of the current integration point. ninpt: The total number of integration points.nnode: The number of nodes per element.
		!ndim: The number of spatial dimensions (typically 2 for 2D elements). !shapef: Array to store the calculated shape functions.
		!dshape: Array to store the derivatives of the shape functions with respect to the natural coordinates.
		
		!! Coordinates of the integration points in the natural coordinate system (r, s)
		data coord24  /1. , 1. ,        &
					  -1. , 1. ,        &
					   1. , -1.,        &
					   -1. , -1./
			
		!Initialize shape functions and their derivatives to zero
		shapef=0.
		dshape=0.
		
		 ! Get the natural coordinates (r, s) for the current integration point
		r=coord24(1,kintk)*gaussCoord
		s=coord24(2,kintk)*gaussCoord
		
		 ! Evaluate the shape functions at the integration point	
		shapef(1)=0.25*(1.+r)*(1.+s)
		shapef(2)=0.25*(1.-r)*(1.+s)
		shapef(3)=0.25*(1.-r)*(1.-s)
		shapef(4)=0.25*(1+r)*(1.-s)
	
		! derivative with respect to r
		dshape(1,1)=0.25*(1.+s)
		dshape(1,2)=-0.25*(1.+s)
		dshape(1,3)=-0.25*(1.-s)
		dshape(1,4)=0.25*(1.-s)	
	
	
	    	! derivative with respect to s
		dshape(2,1)=0.25*(1.+r)
		dshape(2,2)=0.25*(1.-r)
		dshape(2,3)=-0.25*(1.-r)
		dshape(2,4)=-0.25*(1.+r)
	
	
		return
		end
		
!******************************************************
!** Jacobian                                          *
!******************************************************
		subroutine jacobian(ndim,coords,nnode,djac,dshape,xjaci)
		
		include 'aba_param.inc'
		
		real*8 xjac(ndim,ndim),xjaci(ndim,ndim)
		real*8 coords(2,nnode),dshape(ndim,nnode)
 		
 		!The Jacobian matrix transforms the derivatives of the shape functions from the natural coordinate system (r, s) to the global coordinate system (x, y).	
 		!shape functions are defined in the natural coordinates, but the actual calculations (e.g., strains and stresses) are performed in the global coordinates.
	
 		! Initialize the Jacobian matrix and its inverse to zero
		xjac=0.d0
		xjaci=0.d0

		! Calculate the Jacobian matrix	
		do inod= 1, nnode  !The number of nodes per element.
			  do kdim = 1, ndim
				    do jdim = 1, ndim !The number of spatial dimensions (typically 2 for 2D elements).
				    
					      xjac(jdim,kdim) = xjac(jdim,kdim) + dshape(jdim,inod)*coords(kdim,inod) 
					         
				    end do
			  end do 
		end do
		
		! Calculate the determinant of the Jacobian matrix
		djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
		
		! Calculate the inverse(xjaci) of the Jacobian matrix if the determinant(djac) is positive 
		if (djac .gt. 0.) then
			xjaci(1,1)=xjac(2,2)/djac
			xjaci(2,2)=xjac(1,1)/djac
			xjaci(1,2)=-xjac(1,2)/djac
			xjaci(2,1)=-xjac(2,1)/djac			
		end if
		 
		return
		end
!******************************************************
!** B-matrix                                          *
!******************************************************
		subroutine Bmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)
		
		include 'aba_param.inc'
                !implicit none
		
		real*8 xjaci(ndim,ndim),dshape(ndim,4),XYDerivation(ndim,4)
		real*8 bmat(3,2*nnode),bmatT(2*nnode,3)
		
		!Strain-displacement matrix (Size 3*8 for 2D plain stress)
		!The Bmatrix subroutine constructs the B matrix and its transpose. This matrix is used to convert nodal displacements into strains.
		!The B matrix is essential in expressing the strain-displacement relationship within an element. It translates the nodal displacements into strains.
		
		! Initialize XYDerivation and bmat arrays to zero
		XYDerivation=0.d0
		bmat=0.d0
		bmatT=0.d0

		! Compute the derivatives of the shape functions in global coordinates
		!dshape=derivative of shape function
		!do i=1,ndim
			!do j=1,4
				!do k=1,ndim
					XYDerivation=matmul(xjaci,dshape)
				!end do	
			!end do
		!end do
				
		! Construct the B matrix, size=3*8
		bmat(1,1)=XYDerivation(1,1)
		bmat(1,3)=XYDerivation(1,2)
		bmat(1,5)=XYDerivation(1,3)
		bmat(1,7)=XYDerivation(1,4)
		
		bmat(2,2)=XYDerivation(2,1)
		bmat(2,4)=XYDerivation(2,2)
		bmat(2,6)=XYDerivation(2,3)
		bmat(2,8)=XYDerivation(2,4)
		
		bmat(3,1)=XYDerivation(2,1)
		bmat(3,2)=XYDerivation(1,1)
		bmat(3,3)=XYDerivation(2,2)
		bmat(3,4)=XYDerivation(1,2)
		bmat(3,5)=XYDerivation(2,3)
		bmat(3,6)=XYDerivation(1,3)
		bmat(3,7)=XYDerivation(2,4)
		bmat(3,8)=XYDerivation(1,4)
		
		! Construct the transpose of the B matrix
		do i=1,3
			do j=1,2*nnode
				bmatT(j,i)=bmat(i,j)
			end do
		end do
		
		return
		end
		
!*****************************************************************
		  subroutine straininc(ntens,ndof,ndim,nnode,mlvarx,bmat,du,deps)

		  include 'aba_param.inc'
                  !implicit none
		  real*8 deps(3),bmat(3,2*nnode),du(mlvarx,1),xdu(2) 
		  
		  ! Initialize incremental strain vector and temporary displacement increment vector	
		  !du will be provided by abaqus as initial guess and abaqus might be used newton-rapsan mthod to solve the system of linear equaitons
		  deps = 0.d0 !stain increment vector
		  xdu=0.d0  !temporary displacement increment vector
		  
		  !du always provied by abaqus, so abaqus is increaing the displacement i.e. strain step by step until its converge
		  !                        lets abaqus provide du=[du =                    !not given by input file
									!  [0.01],
									!  [0.02],
									!  [0.03],
									!  [0.04],
									!  [0.05],
									!  [0.06],
									!  [0.07],
									!  [0.08]
		  
		  ! Loop over each node in the element	
		  do nodi = 1, nnode                   !number of node per element=4
		  
		           !Compute the starting index for degrees of freedom of the current node   
			   incr_row = (nodi - 1)*ndof  !For each node, it calculates the starting index for the degrees of freedom
			   
			   !Extract displacement increments for the current node
			   do i = 1, ndof         
					xdu(i)= du(i + incr_row,1)
					!here nodof=2
			   		!for nodi=1, incr_row=(1-1)*2=0, first loop
			   		                                 !xdu(1) =  du(1 + incr_row, 1) = du(1, 1) = 0.01  second loop, ndof = 1
			   						 !xdu(2) = du(2 + incr_row, 1) = du(2, 1) = 0.02   second loop, ndof = 2
			   		!for nodi=2, incr_row = (2 - 1) * 2 = 2 first loop			   						  
			   						 !xdu(1) = du(1 + incr_row, 1) = du(3, 1) = 0.03   second loop, ndof = 1
									 !xdu(2) = du(2 + incr_row, 1) = du(4, 1) = 0.04   second loop, ndof = 2
					!for nodi=3, incr_row = (3 - 1) * 2 = 4 first loop
									 !xdu(1) = du(1 + incr_row, 1) = du(5, 1) = 0.05   second loop, ndof = 1
									 !xdu(2) = du(2 + incr_row, 1) = du(6, 1) = 0.06   second loop, ndof = 2
					!for nodi=4, incr_row = (4 - 1) * 2 = 6 first loop
									 !xdu(1) = du(1 + incr_row, 1) = du(7, 1) = 0.07   second loop, ndof = 1
									 !xdu(2) = du(2 + incr_row, 1) = du(8, 1) = 0.08   second loop, ndof = 2			 
			   end do
			   
			   !Compute the contribution of the current node to the incremental strain
	  		   dNidx = bmat(1,1 + incr_row) 
	  		   !this is inside the first do loop
			   !nodi =1, incr_row = (1 - 1) * 2 = 0,   dNidx = bmat(1, 1 + incr_row) = bmat(1, 1)
			   !nodi =2, incr_row = (2 - 1) * 2 = 2,   dNidx = bmat(1, 1 + incr_row) = bmat(1, 3)
			   !nodi =3, incr_row = (3 - 1) * 2 = 4,   dNidx = bmat(1, 1 + incr_row) = bmat(1, 5) 
			   !nodi =4, incr_row = (4 - 1) * 2 = 6,   dNidx = bmat(1, 1 + incr_row) = bmat(1, 7) 
			   
			   dNidy = bmat(2,1 + incr_row + 1) 
			   !this is inside the first do loop
			   !ndof=2
			   !nodi=1, incr_row = (1 - 1) * 2 = 0,   dNidy = bmat(2, 1 + incr_row + 1) = bmat(2, 2)
			   !nodi=2, incr_row = (2 - 1) * 2 = 2,   dNidy = bmat(2, 1 + incr_row + 1) = bmat(2, 4)
			   !nodi=3, incr_row = (3 - 1) * 2 = 4,   dNidy = bmat(2, 1 + incr_row + 1) = bmat(2, 6)   
			   !nodi=4, incr_row = (4 - 1) * 2 = 6,   dNidy = bmat(2, 1 + incr_row + 1) = bmat(2, 8)
				
			   deps(1) = deps(1) + dNidx*xdu(1)                   ! Contribution to ε_xx
			   deps(2) = deps(2) + dNidy*xdu(2)                   ! Contribution to ε_yy
			   deps(3) = deps(3) + dNidy*xdu(1) + dNidx*xdu(2)    ! Contribution to γ_xy
					   
		  end do

		  return
		  end
	
!************************************************************************		
		 subroutine statevar(npt,nsvint,statev,NSVARS,statev_ip,icopy)
	
		  include 'aba_param.inc'
                  !implicit none

		  real*8  statev(NSVARS),statev_ip(nsvint)
		  
		  !The statevar subroutine facilitates the transfer of state variables between the global and local arrays
		  !Before Computation: You need to copy the relevant portion of the global state variables (SVARS) to the local state variables (statevLocal) for the current integration point.
		  !After Computation: You need to copy the updated local state variables (statevLocal) back to the corresponding portion of the global state variables (SVARS).
		  !The subroutine signature includes six parameters, but not all parameters need to be explicitly passed during each call		  
		  !statev: Global array of state variables.
		  !NSVARS: Total number of state variables.
		  !statev_ip: Local array of state variables for the current integration point.
		  !icopy: Copy flag (1 for copying from global to local, 0 for copying from local to global).
		  !kintk maps to npt, SVARS maps to statev, statevLocal maps to statev_ip, 1 or 0 are passed as icopy, other three parameter are explicityly passed.
		  !nsvint: Number of state variables per integration point. nsvint=4, 
		  !npt=1,2,3,4,integration point number or index
		  
		  !isvinc: The calculated starting index for the state variables of the current integration point in the global array statev.
		  !The global state variables array statev stores the state variables for all integration points of all elements sequentially.
		  !Therefore, to access or modify the state variables for a specific integration point, you need to calculate the correct starting index within this global array.
		  !npt - 1): This converts the 1-based integration point index (npt) to a 0-based index, which is common in array indexing.
		  !*nsvint: This scales the 0-based index by the number of state variables per integration point (nsvint) to find the starting position within the global array.
		  
		  isvinc= (npt-1)*nsvint        ! Calculate starting index for the current integration point  npt
                  !isvinc=(1-1)*4=0
                  !isvinc=(2-1)*4=4
                  !isvinc=(3-1)*4=8
                  !isvinc=(4-1)*4=12
                  
		  if (icopy .eq. 1) then
			! Copy from global to local
			do i = 1, nsvint !no of state varaible=4
				  statev_ip(i)=statev(i+isvinc)
				 !thats why variable=16 in input files
				 !it is called inside another do loop in UEL npt=1 to 4, two do loop create 16 elements
				 !local_statev_ip=global-statev
				 !node 1 four variables               !node 2 four variabls                !node 3 four varialbes                !node 4 four variabls
				 !statev_ip(1)=statev(1+0)=statev(1)  !statev_ip(1)=statev(1+4)=statev(5)  !statev_ip(1)=statev(1+8)=statev(9)   !statev_ip(1)=statev(1+12)=statev(13) 
				 !statev_ip(2)=statev(2+0)=statev(2)  !statev_ip(2)=statev(2+4)=statev(6)  !statev_ip(2)=statev(2+8)=statev(10)  !statev_ip(2)=statev(2+12)=statev(14) 
				 !statev_ip(3)=statev(3+0)=statev(3)  !statev_ip(3)=statev(3+4)=statev(7)  !statev_ip(3)=statev(3+8)=statev(11)  !statev_ip(3)=statev(3+12)=statev(15)
				 !statev_ip(4)=statev(4+0)=statev(4)  !statev_ip(4)=statev(4+4)=statev(8)  !statev_ip(4)=statev(4+8)=statev(12)  !statev_ip(4)=statev(4+12)=statev(16) 
			
			end do

		  else       
			! Copy from local to global
			do i = 1, nsvint           !nsvint=4
			
				  statev(i+isvinc)=statev_ip(i)
				  !global-statev= local-statve_ip
				  !node 1 four variables               !node 2 four variabls              !node 3 four varialbes              !node 4 four variabls
				  !statev(1+0)=statev(1)=statve_ip(1)  !statev(1+4)=statev(5)=stave_ip(1) !statev(1+8)=statev(9)=stave_ip(1)  !statev(1+12)=statev(13)=stave_ip(1)
				  !statev(2+0)=statev(2)=statve_ip(2)  !statev(2+4)=statev(6)=stave_ip(2) !statev(2+8)=statev(10)=stave_ip(2) !statev(2+12)=statev(14)=stave_ip(2)
				  !statev(3+0)=statev(3)=statve_ip(3)  !statev(3+4)=statev(7)=stave_ip(3) !statev(3+8)=statev(11)=stave_ip(3) !statev(3+12)=statev(15)=stave_ip(3)
				  !statev(4+0)=statev(4)=statve_ip(4)  !statev(4+4)=statev(8)=stave_ip(4) !statev(4+8)=statev(12)=stave_ip(4) !statev(4+12)=statev(16)=stave_ip(4)
				  
				 
			  
			end do
		  end if

		  return
		  end

!************************************************************************
 subroutine umat(stress,statev,ddsdde,sse,spd,scd,                       &
            rpl,ddsddt,drplde,drpldt,                                   &
            stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,     &
            ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,      &
            celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

    use kvisual
    include 'aba_param.inc'
    !implicit none

    character*80 cmname, filename, text_line, legend, inputfile, odbname, date_time
    real*8  stress(ntens),statev(nstatv),                             &                         
    ddsdde(ntens,ntens),                                                &
    ddsddt(ntens),drplde(ntens),                                        &
    stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),              &
    props(nprops),coords(2),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

	  PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=4,nelem=1,nsdv=4)
	  integer kelem
	  integer NoUEL,unit_number,ierr
	  !character(80) :: text_line
	
	
      Print *, '$$$$$$$$$$$$$$$$UMAT is cought/started$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'  
      ddsdde=0.d0 !no calculation for UMAT
    
      !mapping between UEL and UMAT.
      kelem=int(noel-nelem) 
      !kelem=int(2-1)=1 this is the UEL element number.     
       
      filename = '/home/mmollah/Documents/2024June12/UELtoUMATResultsTransfer.txt'
      legend = 'maping betweeen UEL to UMAT for visualization_10June'
      inputfile = 'SingleElementUEL.inp'
      odbname = 'SingleElementUEL.odb'
      
      unit_number = 300
      open(unit=unit_number, file=filename, status='unknown', action='write', access='sequential', iostat=ierr)
            
      !debugging step whether the txt file is created or not
      !If there is an error opening the file, an error message is printed and the subroutine stops execution.
      if (ierr /= 0) then
          print *, 'Error opening file: ', filename
          stop
      endif
      
      write(unit_number, '(A)') trim(legend)
      write(unit_number, '(A, A, A)') 'Input File: ', trim(inputfile)
      write(unit_number, '(A, A, A)') 'ODB Name: ', trim(odbname)
      !write(unit_number, '(A, A, A)') 'Date and Time: ', trim(date_time)
      write(unit_number,'(A)') 'kelem=UEL specified element number, kinpt=npt=specific integration point number,'
      write(unit_number, '(A)') 'Note: K1=1 means S11=SDVI1, K1=2 means S22=SDVI2, K1=3 means S33, K1=4 means S12'
      write(unit_number, '(A)') '---------------------------------------------'
      write(unit_number, '(A)') '    kelem,             npt,                S11            S22            S33            S12'
      
      !This loop iterates over each integration point and each state variable, copying the values from sigout to statev.
      do kintk=1, npt 
	      do k1 = 1, nsdv
		statev(k1) = sigout(kelem,kintk,k1)
		!write(unit_number, '(I8, 13X, I4, 8X, 4(1X, E15.6))') kelem, npt, statev(1), statev(2), statev(3), statev(4)
	      end do
	      write(unit_number, '(I8, 13X, I4, 8X, 4(1X, E15.6))') kelem, kintk, statev(1), statev(2), statev(3), statev(4)
      end do
      
      ! Write "Hello World" to the file
      write(unit_number, '(A)') 'Hello World_UMAT-June10 nsdv'
      Print *, 'hellow world from UMAT'
      !Close the file
      close(unit_number)
      Print *, '$$$$$$$$$$$$$$$$UMAT is ended$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      return
      end 
