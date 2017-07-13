/*
#   Apply constraints on tapers using simple derivatives
#
#     Copyright (C) 2013-2017  Andrew J. Barbour *
#
#     Robert L. Parker authored the original matlab algorithm
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
#
*/

#include <stdio.h>
#include <R.h>  
#include <Rinternals.h>
#include <Rdefines.h>

// c implementation of original RLP constraint filter
SEXP rlp_constrain(SEXP R_ntaps, SEXP R_maxslope)
{
    // copy semantics was leading to changes in env vars
    // so we now protect and duplicate
    // http://adv-r.had.co.nz/C-interface.html
    SEXP R_ntaps_copy = PROTECT(duplicate(R_ntaps));
    double * c_ntaps=REAL(R_ntaps_copy);
    UNPROTECT(1);
    int maxslope=REAL(R_maxslope)[0];
    int ssize=LENGTH(R_ntaps_copy);
    SEXP ntap_con;
    int state, i, im, msize = ssize - 1;
    double slope;

    if (maxslope <= 0)
        Rf_error( "max slope must greater than zero" );

    // APPLY CONSTRAINTS
    //     FORWARD:
    state = 0;
    for (i = 1; i < ssize; i++){
        im = i - 1;
        if (state == 0){
          slope = c_ntaps[ i ] - c_ntaps[ im ];
    			if (slope >= maxslope) {
    				state = 1;
    				c_ntaps[ i ] = c_ntaps[ im ] + maxslope;
    			}
        } else {
    			if (c_ntaps[ i ] >= c_ntaps[ im ] + maxslope) {
    				c_ntaps[ i ] = c_ntaps[ im ] + maxslope;
    			} else {
    				state = 0;
    			}
        }
    }
    // and BACKWARD:
    state = 0;
    for (i = msize; i >= 1; i--){
        im = i - 1;
        if (state == 0){
            slope = c_ntaps[ im ] - c_ntaps[ i ];
            if (slope >= maxslope) {
            	state = 1;
            	c_ntaps[ im ] = c_ntaps[ i ] + maxslope;
            }
        } else {
            if (c_ntaps[ im ] >= c_ntaps[ i ] + maxslope) {
            	c_ntaps[ im ] = c_ntaps[ i ] + maxslope;
            } else {
            	state = 0;
            }
        }
    }
    PROTECT(ntap_con = allocVector(REALSXP, ssize));
    for (i = 0; i < ssize; i++){
        REAL(ntap_con)[i] = c_ntaps[i];
    }
    UNPROTECT(1);
    return(ntap_con);
}

