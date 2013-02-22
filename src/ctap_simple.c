/*
#   psd: 
#
#   Copyright (C) 2013  Andrew J. Barbour andy.barbour @ gmail.com
#
#   Robert L. Parker authored the original algorithm.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
// R specific libs:
#include <R.h>  
#include <Rinternals.h>
#include <Rdefines.h>

/* Sending R integer vectors to C using .Call */
SEXP rlp_constrain_tapers(SEXP R_ntaps, SEXP R_maxslope)
{
    // http://r.789695.n4.nabble.com/protect-unprotect-howto-in-C-code-td911033.html
    // http://stackoverflow.com/questions/4106174/where-can-i-learn-to-how-to-write-c-code-to-speed-up-slow-r-functions
    // PROTECTION [ ]
    //
    //PROTECT(R_maxslope = AS_INTEGER(R_maxslope));
    //PROTECT(c_ntaps = AS_NUMERIC(R_ntaps));
    double * c_ntaps=REAL(R_ntaps);
    int maxslope=REAL(R_maxslope)[0];
    int ssize=LENGTH(R_ntaps);
    SEXP ntap_con;
    //--------------------------------------------//
    //
    //--------------------------------------------//
    int state, i, l;
    // forward indices
    int i_min_f = 1, i_max_f = ssize - 1;
    // reverse indices
    int i_min_r = i_max_f - 1, i_max_r = 0;
    // num-tapers and slopes
    double slope;
    int ntap_c, ntap_p, ntap_fix;
    if(maxslope <= 0)
        Rf_error( "max slope must greater than zero" );
    //
    double working_space[ssize];
    //
    // RLPs algorithm
    //  #  Scan forward to bound slopes >= 1
    // state<-0
    // for ( j  in  2:nf ) {
    //   slope <- slopes[j]
    //   if (state == 0) {
    //     if (slope >= 1 ) {
    //       state <- 1
    //       kopt[j] <- kopt[j-1]+1
    //     }
    //   } else {
    //     if (kopt[j] >= kopt[j-1]+1) {
    //       kopt[j] <- kopt[j-1]+1
    //     } else {
    //       state <- 0
    //     }
    //   }
    // }
    //  #  Scan backward to bound slopes >= -1
    // 	state <- 0
    // 	for ( j  in  nf:2 ) {
    //         if (state == 0) {
    // 			slope <- kopt[j-1]-kopt[j]
    // 			if (slope >= 1) {
    // 				state <- 1
    // 				kopt[j-1] <- kopt[j]+1
    // 			}
    //         } else {
    // 			if (kopt[j-1] >= kopt[j]+1) {
    // 				kopt[j-1] <- kopt[j]+1
    // 			} else {
    // 				state <- 0
    // 			}
    //         }
    // 	}
    // APPLY CONSTRAINTS
    //     FORWARD:
    state = 0;
    for (i = i_min_f; i <= i_max_f; i++){
        if (state == 0){
            slope = c_ntaps[i] - c_ntaps[i-1];
			if (slope >= maxslope) {
				state = 1;
				//printf("sub f-a\n");
				c_ntaps[i] = c_ntaps[i-1] + maxslope; // was orig 1
			}
        } else {
			if (c_ntaps[i] >= c_ntaps[i-1] + maxslope) {
			    //printf("sub f-b\n");
				c_ntaps[i] = c_ntaps[i-1] + maxslope;
			} else {
				state = 0;
			}
        }
	    //printf("f %i %f %i %f\n", i, slope, state, c_ntaps[i]);
    }
    // and BACKWARD:
    for (i = i_min_r; i >= i_max_r; i--){
        if (state == 0){
            slope = c_ntaps[i-1] - c_ntaps[i];
			if (slope >= maxslope) {
				state = 1;
				//printf("sub r-a\n");
				c_ntaps[i-1] = c_ntaps[i] + maxslope; // was orig 1
			}
        } else {
			if (c_ntaps[i-1] >= c_ntaps[i] + maxslope) {
			    //printf("sub r-b\n");
				c_ntaps[i-1] = c_ntaps[i] + maxslope;
			} else {
				state = 0;
			}
        }
	    //printf("f %i %f %i %f\n", i, slope, state, c_ntaps[i]);
    }
    // use t
    PROTECT(ntap_con = allocVector(REALSXP, ssize));
    // C arrays begin index zero
    for (i = 0; i < ssize; i++){
        REAL(ntap_con)[i] = c_ntaps[i];
    }
    UNPROTECT(1);
    return(ntap_con);
}

// http://www.sfu.ca/~sblay/R-C-interface.txt
// V.:
/* useCall3.c                                    */
/* Getting an integer vector from C using .Call  */
SEXP setInt() {
    SEXP myint;
    int *p_myint; 
    int len = 5;
    PROTECT(myint = NEW_INTEGER(len));  // Allocating storage space
    p_myint = INTEGER_POINTER(myint);
    p_myint[0] = 7;
    UNPROTECT(1);
    return myint;
}
// to work with real numbers, replace int with double and INTEGER with NUMERIC
