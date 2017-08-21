#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwo2xsf_old.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# set locales to C
LANG=C 
LC_ALL=C
export LANG LC_ALL

#
# pwo2xsf_old.sh: PW-output--to--XSF conversion; 
#                 Used for PWscf versions < 1.3
#
# Usage: pwo2xsf_old.sh [options] pw-output-file
#
# Written by Tone Kokalj on Tue May  8 20:43:44 2001
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cat > pwo2xsfUsage.$$ <<EOF

 Usage: pwo2xsf.sh [options] [pw-output-file]

 Options are:

               pwo2xsf.sh --initcoor|-ic [pw-output-file]
                             Will extract the initial (i.e. input) 
                             ionic coordinates.

               pwo2xsf.sh --latestcoor|-lc [pw-output-file]
                             Will extract latest estimation of ionic 
                             coordinates from pw-output file. The coordinates 
                             can be either taken from "Search of equilibrium 
                             positions" record or from "Final estimate of 
                             positions" record.

               pwo2xsf.sh --optcoor|-oc [pw-output-file]
			     Similar to "--latestcoor", but extract just the 
			     optimized coordinates from "Final estimate of 
                             positions" record.

	       pwo2xsf.sh --animxsf|-a [pw-output-file1] ...
			     Similar to "--latestcoor", but extract the
			     coordinates from all ionic steps and make
			     a AXSF file for animation.

               pwo2xsf.sh -r option [pw-output-file]
			     one must specify i.e. ityp->nat conversion, 
			     and the corresponding data are written to file 
			     nuclei.charges. The -r flag deletes this file.
                             Here "option" is one of the above options
			     (i.e -lc|-oc|-a).
EOF

pwoCleanFiles() {
    for file in pwo2xsfUsage.$$ pw.$$ pwo2xsf.xsf1 pwo2xsf.xsf2
    do
      if test -f $file; then rm -f $file; fi
    done
}
pwoExit() {
    # Usage: $0 status
    pwoCleanFiles
    exit $status
}
pwoUsage() {
    if [ $1 ]; then
	echo "
Usage: $2
"
	pwoExit 1
    fi
}

#######################################
if [ "$XCRYSDEN_TOPDIR" != "" ]; then
    . $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh
else
    echo "
    ERROR: cannot convert to XSF, because loading of pwLib_old.sh failed
"
    exit 1
fi
#######################################

#############
BOHR=0.529177
#############

# ---------------------------------------------------------------------------
# read PW-output file and print the XSF file according to specified flags
pwoPrintCoor() {
    #set -x
    pwoUsage "$# -lt 1" \
	"$0 --latestcoor|-lc [pw-output-file]   or   $0 --optcoor|-oc [pw-output-file]"
   
    case $1 in
	--latestcoor|-lc) type=lc; shift;;
	--optcoor|-oc)    type=oc; shift;;
    esac

    if [ $# -eq 0 ]; then
	pwNucleiCharges -  pw.$$
	cat - >> pw.$$
    else
	pwNucleiCharges $1  pw.$$
	for i in `ForLoop 1 $#`
	do
	    cat $1 >> pw.$$
    	    shift	    
	done
    fi
    inp=pw.$$

    cat "$inp" | awk -v bohr=$BOHR -v t=$type -- '
function isInt(__num) {
   if ( __num ~ /^[0-9]+$/ ) return 1;
   else return 0;
} 
function PrintItyp(__ityp) {
   if ( __ityp ~ /^[0-9]+$/ ) return atn[ __ityp ];
   else return __ityp;
}
BEGIN { 
    line=0; optc=0; 
    for (i=0; i<100; i++) atn[i]=i;
    getline;
    ntyp=$1;
    for (i=0; i<ntyp; i++) {
	getline;
	atn[$1] = $2;
    }
}

/bravais-lattice index/       { ibrav=$NF; }
/celldm\(1\)=/                { a0=$2*bohr; scale=a0; }
/crystal axes:/               {
    for (i=0; i<3; i++) {
	getline;
	for (j=4; j<7; j++) {
	    vec[i,j-4] = $j * a0;
	}
    }
    printf " DIM-GROUP\n 3 1\n PRIMVEC\n";
    printf "%15.10f %15.10f %15.10f\n", vec[0,0], vec[0,1], vec[0,2];
    printf "%15.10f %15.10f %15.10f\n", vec[1,0], vec[1,1], vec[1,2];
    printf "%15.10f %15.10f %15.10f\n", vec[2,0], vec[2,1], vec[2,2];
}
/number of atoms/             { 
    nat=$NF; 
}
/Final estimate of positions/ { 
    optc=1; line=1; 	    
    printf " PRIMCOORD\n %d 1\n", nat;
    next; 
}
/Carthesian coordinates/ {
    if (optc==1 && line==1) {
	next;
    }
}
/Search of equilibrium positions/ { if ( t == "lc" ) {
    getline; getline; 
    if ($1 ==  "ATOMIC_POSITIONS" ) {
       if ( $2 ~ /alat/ ) scale=a0;
       else if ( $2 ~ /angstrom/ ) scale=1.0;
       else if ( $2 ~ /bohr/ ) scale=bohr;
       else if ( $2 ~ /crystal/ ) {
          print "ERROR: pwo2xsf.sh does not work for \"crystal\" coordinates";
          exit 1;
       }
       getline;
    }
    # maybe current line is: "Carthesian coordinates"
    if ( NF != 4 ) getline;
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=scale*$1; y[i]=scale*$2; z[i]=scale*$3; ityp[i]=$4;
        } else { 
           x[i]=scale*$2; y[i]=scale*$3; z[i]=scale*$4; ityp[i]=$1;
        }
	getline;
    }
  }
}
/Entering Dynamics;/ {
    getline; getline; getline; nstep++;
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=a0*$1; y[i]=a0*$2; z[i]=a0*$3; ityp[i]=$4;
        } else { 
           x[i]=a0*$2; y[i]=a0*$3; z[i]=a0*$4; ityp[i]=$1;
        }
	getline;
    }
}

/a*/ { 
    if (line == 1) {
       if ($1 ==  "ATOMIC_POSITIONS" ) {
          if ( $2 ~ /alat/ ) scale=a0;
          else if ( $2 ~ /angstrom/ ) scale=1.0;
          else if ( $2 ~ /bohr/ ) scale=bohr;
          else if ( $2 ~ /crystal/ ) {
             print "ERROR: pwo2xsf.sh does not work for \"crystal\" coordinates";
             exit 1;
          }
          getline;
       }

       if (NF != 4) exit 0;	    	    

       if (isInt($4)) printf "% 3d   % 15.10f  % 15.10f  % 15.10f\n",
                             atn[$4], scale*$1, scale*$2, scale*$3; 
       else printf "% 3s   % 15.10f  % 15.10f  % 15.10f\n",
                   $1, scale*$2, scale*$3, scale*$4;
    } 
}
END { if ( t == "lc" && line == 0) {
    printf " PRIMCOORD\n %d 1\n", nat;
    for(i=1; i<=nat; i++)
       printf "% 3s   % 15.10f  % 15.10f  % 15.10f\n", 
	    PrintItyp(ityp[i]), x[i], y[i], z[i];
       }
    }' > pwo2xsf.xsf_out
    cat pwo2xsf.xsf_out
}

# ---------------------------------------------------------------------------
# read PW-output file and print the animated-XSF file
pwoAnimXSF() {
    #set -x
    pwoUsage "$# -lt 1" "$0 --animxsf|-a [pw-output-file1] ..."
    
    only_init=0
    case $1 in
	--inicoor|-ic) only_init=1;;
    esac

    if [ $# -eq 1 ]; then
	pwNucleiCharges -  pw.$$    
	cat - >> pw.$$
    else
	pwNucleiCharges $2  pw.$$
	# if "pw-output-file" is glob expresion -> must merge all outputs
	for i in `ForLoop 2 $#`
	do
	    cat $2 >> pw.$$
	    shift	    
	done
    fi
    inp=pw.$$

    nstep=`egrep "Final estimate of positions|Search of equilibrium positions|Entering Dynamics;" $inp | wc | awk '{print $1}'`
    # add also initial coordinates
    nstep=`expr $nstep + 1`

    cat "$inp" | awk -v bohr=$BOHR -v astep=$nstep -v onlyinit=$only_init -- '
function isInt(__num) {
   if ( __num ~ /^[0-9]+$/ ) return 1;
   else return 0;
}
function PrintItyp(__ityp) {
   if ( __ityp ~ /^[0-9]+$/ ) return atn[ __ityp ];
   else return __ityp;
}

function PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz) {
    if (onlyinit) {
       print " PRIMCOORD";
    } else {
       print " PRIMCOORD", nstep;
    }
    print nat, 1;
    for(i=1; i<=nat; i++) {
	printf "% 3s ", PrintItyp(ityp[i]); 
        printf "% 15.10f  % 15.10f  % 15.10f   % 15.10f  % 15.10f  % 15.10f\n", x[i], y[i], z[i], fx[i], fy[i], fz[i];
    }
}
function GetInitCoor(nstep, nat, scale, ityp, x, y, z) {
    for(i=1; i<=nat; i++) {
        ityp[i]=$2;
	split($0,rec,"("); split(rec[3],coor," ");
	x[i]= scale*coor[1]; y[i]=scale*coor[2]; z[i]=scale*coor[3];
	getline;
    }
}
function GetPrimCoor(nstep, nat, scale, ityp, x, y, z) {
    for(i=1; i<=nat; i++) {
        if (isInt($4)) {
           x[i]=scale*$1; y[i]=scale*$2; z[i]=scale*$3; ityp[i]=$4;
        } else { 
           x[i]=scale*$2; y[i]=scale*$3; z[i]=scale*$4; ityp[i]=$1;
        }
	getline;
    }
}
function GetForces(nstep, nat, ityp, fx, fy, fz) {
    for(i=1; i<=nat; i++) {
	#if (nstep==1) ityp[i]=$4;
	fx[i]=$7; fy[i]=$8; fz[i]=$9;
	getline;
    }
}

BEGIN { 
    nstep=1; print_celldm=1; 
    for (i=0; i<100; i++) atn[i]=i;
    getline;
    ntyp=$1;
    for (i=0; i<ntyp; i++) {
	getline;
	atn[$1] = $2;
    }
}

/bravais-lattice index/       { ibrav=$NF; }

/celldm\(1\)=/                { a0=$2*bohr; scale=a0; }

/crystal axes:/               {
    if (print_celldm) {
	for (i=0; i<3; i++) {
	    getline;
	    for (j=4; j<7; j++) {
		vec[i,j-4] = $j * a0;
	    }
	}
        #---
        # this is now donw posteriori (see below)
        #if (!onlyinit) printf " ANIMSTEPS %d\n", astep
        #---
	printf " DIM-GROUP\n 3 1\n PRIMVEC\n";
	printf "%15.10f %15.10f %15.10f\n", vec[0,0], vec[0,1], vec[0,2];
	printf "%15.10f %15.10f %15.10f\n", vec[1,0], vec[1,1], vec[1,2];
	printf "%15.10f %15.10f %15.10f\n", vec[2,0], vec[2,1], vec[2,2];
	print_celldm=0;
    }
}

/number of atoms/  { nat=$NF; }

/Carthesian axes|Cartesian axes/  { 
    getline; getline; getline;
    if (nstep == 1) GetInitCoor(nstep, nat, a0, ityp, x, y, z);
    if (onlyinit == 1) {
       PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz);
       exit;
    }
}
/Final estimate of positions/     { 
    getline; nstep++;

    if ($1 ==  "ATOMIC_POSITIONS" ) {
       if ( $2 ~ /alat/ ) scale=a0;
       else if ( $2 ~ /angstrom/ ) scale=1.0;
       else if ( $2 ~ /bohr/ ) scale=bohr;
       else if ( $2 ~ /crystal/ ) {
          print "ERROR: pwo2xsf.sh does not work for \"crystal\" coordinates";
          exit 1;
       }
       getline;
    }

    # maybe current line is: "Carthesian coordinates"
    if ( NF != 4 ) getline;
    GetPrimCoor(nstep, nat, scale, ityp, x, y, z);
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz); 
}

/Search of equilibrium positions/ { 
    getline; getline; nstep++;

    if ($1 ==  "ATOMIC_POSITIONS" ) {
       if ( $2 ~ /alat/ ) scale=a0;
       else if ( $2 ~ /angstrom/ ) scale=1.0;
       else if ( $2 ~ /bohr/ ) scale=bohr;
       else if ( $2 ~ /crystal/ ) {
          print "ERROR: pwo2xsf.sh does not work for \"crystal\" coordinates";
          exit 1;
       }
       getline;
    }


    # maybe current line is: "Carthesian coordinates"
    if ( NF != 4 ) getline;
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    GetPrimCoor(nstep, nat, scale, ityp, x, y, z);
}
/Entering Dynamics;/ {
    getline; getline; getline; nstep++;
    for (i=1; i<=nat; i++) {
	fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;
    }
    GetPrimCoor(nstep, nat, a0, ityp, x, y, z);
}
/Forces acting on atoms/          { 
    getline; getline; 
    GetForces(nstep, nat, ityp, fx, fy, fz); 
    PrintPrimCoor(nstep, nat, atn, ityp, x, y, z, fx, fy, fz); 
}' > pwo2xsf.xsf_out

    if [ $only_init -eq 0 ]; then
    # Assign the number of ANIMSTEPS here. The reason is that the
    # output file (queue runs) is the result of several job runs, then
    # some of them might be terminated on the "wrong" place, and the
    # initial ANIMSTEPS might be wrong. The most secure way is to extract the 
    # sequential digit from the last "PRIMCOORD id" record.
	nsteps=`grep PRIMCOORD pwo2xsf.xsf_out | tail -1 | awk '{print $2}'`    
	echo " ANIMSTEPS $nsteps" > pwo2xsf.xsf1
	cp pwo2xsf.xsf_out pwo2xsf.xsf2
	cat pwo2xsf.xsf1 pwo2xsf.xsf2 > pwo2xsf.xsf_out
    fi
    
    cat pwo2xsf.xsf_out
}


#######################################################################
####                              MAIN                              ###
#######################################################################
if [ $# -eq 0 ]; then
    cat pwo2xsfUsage.$$
    pwoExit 1
fi

r=0
if [ "$1" = "-r" ]; then
    r=1
    shift
fi

case $1 in
    --inicoor|-ic)    pwoAnimXSF $@;;
    --latestcoor|-lc) pwoPrintCoor $@;;
    --optcoor|-oc)    pwoPrintCoor $@;;
    --animxsf|-a)     pwoAnimXSF $@;;    
    *) cat pwo2xsfUsage.$$; pwoExit;;
esac

if [ $r -eq 1 ]; then
    rm nuclei.charges
fi

pwoExit 0