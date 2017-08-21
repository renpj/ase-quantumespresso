#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwo2xsf_neb.awk
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# Purpose: extract all coordinates and forces from the NEB (path) PW.out file for PWscf-v2.0 or latter
# BEWARE:  NEB (path) and variable-cell is currently not supported
#
# Written by Tone Kokalj on Thu Jan  5 14:16:36 CET 2006
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

function PrintPrimCoor(istep, nat, atom, x, y, z, fx, fy, fz) {
  print " PRIMCOORD", istep;  
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


function make_error(message,status) {
  printf "ERROR: %s\n", message;
  error_status=status;
  exit status;
}


BEGIN {
  bohr=0.529177;
  istep=1;
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
    #for (j=4; j<7; j++) v[i,j-4] = $j * a0;      
  }
  if (istep==1) {
    PrintSpecies();
    PrintPrimVec(is_vc,istep,v);
  }
}


/Cartesian axes/  { 
  # read INITIAL coordinates
  getline; getline; getline;
  GetInitCoor(nat, a0, atom, x, y, z);
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
  
  PrintPrimCoor(istep, nat, atom, x, y, z, fx, fy, fz); 
  istep++;
}

