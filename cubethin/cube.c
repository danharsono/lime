/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool 
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2013, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"

/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->radius		= 2000*AU;
  par->minScale	   	= 0.1*AU;
  par->pIntensity    	= 60000;
  par->sinkPoints    	= 6000;
  par->dust		= "jena_thin_e6.tab";
  par->moldatfile[0] 	= "hco+@xpol.dat";
  par->antialias	= 4;
  par->sampling		= 1;

  par->outputfile 	= "populations.pop";
  par->binoutputfile 	= "restart.pop";
  par->gridfile		= "grid.vtk";

/* 
 * Definitions for image #0. Add blocks for additional images.
 */
//  img[0].nchan			= 60;		  // Number of channels
//  img[0].velres			= 500.;       // Channel resolution in m/s
//  img[0].trans			= 3;          // zero-indexed J quantum number
  img[0].freq           =   115e9;
  img[0].pxls		= 512;	      // Pixels per dimension
  img[0].imgres		= 0.1;		  // Resolution in arc seconds
  img[0].theta		= 0.0;		  // 0: face-on, pi/2: edge-on
  img[0].distance	= 100*PC;	  // source distance in m
  img[0].source_vel	= 0;          // source velocity in m/s
  img[0].unit		= 1;		  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename	= "cube.fits";	// Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){	
    /*
    * Define variable for radial coordinate
    */
    double r;
    /* 
    * Calculate radial distance from origin
    */
    if(x<500*AU && x>-500*AU && y<500*AU && y>-500*AU && z<500*AU && z>-500*AU){
        density[0] = 1e6;
    } else {
    density[0] = 0.;
    }
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
    int i,x0=0;
    double r;
    /* 
    * Array containing temperatures as a function of radial 
    * distance from origin (this is an example of a tabulated model)
    */
    if(x<500*AU && x>-500*AU && y<500*AU && y>-500*AU && z<500*AU && z>-500*AU){
        temperature[0] = 60;
    } else {
        temperature[0] = 2.7;
}
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
/* 
 * Here we use a constant abundance. Could be a 
 * function of (x,y,z).
 */
  abundance[0] = 1.e-9;
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
/* 
 * 200 m/s as the doppler b-parameter. This
 * can be a function of (x,y,z) as well.
 * Note that *doppler is a pointer, not an array. 
 * Remember the * in front of doppler.
 */
  *doppler = 200.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
/*
 * Variables for spherical coordinates
 */
  double R, phi,r,theta;
/*
 * Transform Cartesian coordinates into spherical coordinates
 */
  R=sqrt(x*x+y*y+z*z);
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
/*
 * Free-fall velocity in the radial direction onto a central 
 * mass of 1.0 solar mass
 */  
  r=-sqrt(2*6.67e-11*1.989e30/R);
/*
 * Vector transformation back into Cartesian basis
 */
  vel[0]=r*sin(theta)*cos(phi);
  vel[1]=r*sin(theta)*sin(phi);
  vel[2]=r*cos(theta);
}

/******************************************************************************/


