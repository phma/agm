/******************************************************/
/*                                                    */
/* ldecimal.h - lossless decimal representation       */
/*                                                    */
/******************************************************/
/* Copyright 2015,2017,2020 Pierre Abbat.
 * This file is part of AGM.
 */

#include <string>

std::string ldecimal(double x,double toler=0,bool noexp=false);
/* Returns the shortest decimal representation necessary for
 * the double read back in to be equal to the double written.
 * If toler>0, returns the shortest representation of a number
 * that is within toler of x.
 */
