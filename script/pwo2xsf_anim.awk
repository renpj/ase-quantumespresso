#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwo2xsf_anim.awk
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# Purpose: extract either INITIAL or all coordinates for PWscf-v2.0 or latter
#
# Written by Tone Kokalj on Mon Feb  9 12:48:10 CET 2004
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function PrintSpecies() {
  if (dimen==0)      print "MOLECULE";
  else if (dimen==1) print "POLYMER";
  else if (dimen==2) print "SLAB";
  else               print "CRYSTAL";
}		  

function PrintPrimVec(is_vc,ith,vec) {
  if (!is_vc) printf "PRIMVEC\n";
  else        printf "PRIMVEC %d\n",ith;
  printf "  %15.10f %15.10f %15.10f\n", v[0,0], v[0,1], v[0,2];
  printf "  %15.10f %15.10f %15.10f\n", v[1,0], v[1,1], v[1,2];
  printf "  %15.10f %15.10f %15.10f\n", v[2,0], v[2,1], v[2,2];  
}

function PrintPrimCoor(onlyinit,istep, nat, atom, x, y, z, fx, fy, fz) {
  if (onlyinit) {
    print " PRIMCOORD";
  } else {
    print " PRIMCOORD", istep;
  }
  print nat, 1;
  
  for(i=0; i<nat; i++) {
    printf "  %3s    % 15.10f  % 15.10f  % 15.10f    % 15.10f  % 15.10f  % 15.10f\n", atom[i], x[i], y[i], z[i], fx[i], fy[i], fz[i];
  }
}

function GetInitCoor(nat, scale, atom, x, y, z) {
  for(i=0; i<nat; i++) {
    atom[i]=$2;
    split($0,rec,"("); split(rec[3],coor," ");
    x[i]= scale*coor[1]; y[i]=scale*coor[2]; z[i]=scale*coor[3];
    getline;
  }
}

function CrysToCartCoor(i,v,a,b,c) {
  # Crystal --> Cartesian (ANGSTROM units) conversion
  x[i] = v[0,0]*a + v[1,0]*b + v[2,0]*c;
  y[i] = v[0,1]*a + v[1,1]*b + v[2,1]*c;
  z[i] = v[0,2]*a + v[1,2]*b + v[2,2]*c;
}

function make_error(message,status) {
  printf "ERROR: %s\n", message;
  error_status=status;
  exit status;
}


BEGIN {
  bohr=0.529177;
  istep=1;
  error_status=0;
  if (nvec>1 || (nvec==1 && ncoor==2)) {
    is_vc=1; # variable-cell
  } else {
    is_vc=0;
  }
}


/celldm\(1\)=/    { a0=$2*bohr; scale=a0; l_scale=a0; }
/number of atoms/ { nat=$NF; }


/crystal axes:/   {
  # read the lattice-vectors
  for (i=0; i<3; i++) {
    getline;
    split($0,rec,/\(/);
    split(rec[3],vecstr," ");
    for (j=1; j<4; j++) v[i,j-1] = vecstr[j] * a0;    
  }
  if (istep==1) {
    PrintSpecies();
    PrintPrimVec(is_vc,istep,v);
  }
}


/Cartesian axes/  { 
  # read INITIAL coordinates
  getline; getline; getline;
  if (istep == 1) GetInitCoor(nat, a0, atom, x, y, z);
}


$1 == "CELL_PARAMETERS" {
  # read the lattice-vectors (type=LATEST and OPTIMIZED)
  ff=l_scale;
  if      ( $2 ~ /alat/ )     ff=a0;
  else if ( $2 ~ /angstrom/ ) ff=1.0;
  else if ( $2 ~ /bohr/ )     ff=bohr;
  for (i=0; i<3; i++) {
    getline;
    if (NF != 3) make_error("error reading CELL_PARAMETERS records",1);
    for (j=1; j<4; j++) {
      v[i,j-1] = ff * $j;       
    }
  }
  if (is_vc) PrintPrimVec(is_vc,istep,v);
}  


$1 == "ATOMIC_POSITIONS" {
  # read atomic positions
  crystal_coor=0;
  if      ( $2 ~ /alat/ )     scale=a0;
  else if ( $2 ~ /angstrom/ ) scale=1.0;
  else if ( $2 ~ /bohr/ )     scale=bohr;
  else if ( $2 ~ /crystal/ ) {
    scale=1.0;
    crystal_coor=1;
  }
  for(i=0; i<nat; i++) {
    getline;
    if (NF != 4 && NF != 7) make_error("error reading ATOMIC_POSITIONS records",1);
    atom[i]=$1; 
    a=scale*$2;
    b=scale*$3;
    c=scale*$4; 
    if (crystal_coor) {
      CrysToCartCoor(i,v,a,b,c); 
    } else {
      x[i]=a; y[i]=b; z[i]=c; 
    }
  }
}


/Forces acting on atoms/ { 
  # read forces
  while ( $1 != "atom" ) {
    if (getline <= 0) {
      # unexpected EOF or error; set forces to zero
      for(i=0; i<nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
      }
      PrintPrimCoor(onlyinit,istep, nat, atom, x, y, z, fx, fy, fz); 
      next;
    }
  }

  for(i=0; i<nat; i++) {
    if (NF != 9) make_error("error reading Forces records",1);
    fx[i]=$7; fy[i]=$8; fz[i]=$9;
    getline;
  }

  if (onlyinit) exit 0;
  else PrintPrimCoor(onlyinit,istep, nat, atom, x, y, z, fx, fy, fz); 
  istep++;
}


END {
  if (error_status == 0 && onlyinit == 1) {
    PrintPrimCoor(onlyinit,istep, nat, atom, x, y, z, fx, fy, fz);
  }
}
