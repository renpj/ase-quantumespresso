#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwo2xsf_opt.awk
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# Purpose: extract the LATEST or OPTIMIZED coordinates for PWscf-v2.0 or latter
#
# Written by Tone Kokalj on Mon Feb  9 12:48:10 CET 2004
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function PrintSpecies() {
  if (dimen==0)      print "MOLECULE";
  else if (dimen==1) print "POLYMER";
  else if (dimen==2) print "SLAB";
  else               print "CRYSTAL";
}		  
		  
function CheckAtoms() {
  if (nat < 1) {
    print "ERROR: no atoms found";
    error_status=1;
    exit 1;
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
  nat=0; 
  opt_coor_found=0; 
  error_status=0;
  bohr=0.529177
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
}


$1 == "CELL_PARAMETERS" {
  # read the lattice-vectors (type=OPTIMIZED)
  opt_coor_found=1;
  ff=l_scale;
  if ( $2 ~ /alat/ )          ff=a0;
  else if ( $2 ~ /angstrom/ ) ff=1.0;
  else if ( $2 ~ /bohr/ )     ff=bohr;
  CheckAtoms();
  for (i=0; i<3; i++) {
    getline;
    if (NF != 3) make_error("error reading CELL_PARAMETERS records",1);
    for (j=1; j<4; j++) {
      v[i,j-1] = ff * $j;       
    }
  }
}    


$1 == "ATOMIC_POSITIONS" {
  crystal_coor=0;
  if ( $2 ~ /alat/ )          scale=a0;
  else if ( $2 ~ /angstrom/ ) scale=1.0;
  else if ( $2 ~ /bohr/ )     scale=bohr;
  else if ( $2 ~ /crystal/ ) {
    scale=1.0;
    crystal_coor=1;
  }
  CheckAtoms();
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
      next;
    }
  }

  for(i=0; i<nat; i++) {
    if (NF != 9) make_error("error reading Forces records",1);
    fx[i]=$7; fy[i]=$8; fz[i]=$9;
    getline;
  }    
}

/Begin final coordinates/ {
    opt_coor_found=1;
}
/Final estimate of positions/ {
  # the PWscf v.1.3.0 has still this record
  opt_coor_found=1;
}


END {
  if (error_status==0) {
    if (t=="OPTIMIZED" && !opt_coor_found) {
      print "ERROR: no optimized coordinates";
      exit 1;
    }      
    PrintSpecies();
    printf "PRIMVEC\n";
    printf "  %15.10f %15.10f %15.10f\n", v[0,0], v[0,1], v[0,2];
    printf "  %15.10f %15.10f %15.10f\n", v[1,0], v[1,1], v[1,2];
    printf "  %15.10f %15.10f %15.10f\n", v[2,0], v[2,1], v[2,2];
    printf "PRIMCOORD\n %d 1\n", nat;

    for(i=0; i<nat; i++)
      printf "%  3s   % 15.10f  % 15.10f  % 15.10f   % 15.10f  % 15.10f  % 15.10f\n", atom[i], x[i], y[i], z[i], fx[i], fy[i], fz[i];
  }
}
