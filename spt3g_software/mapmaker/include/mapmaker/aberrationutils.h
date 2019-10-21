#pragma once
#include <G3Map.h>
#include <G3Vector.h>
#include <coordinateutils/G3SkyMap.h>
void aberrate_delta_alpha(G3MapVectorDouble alpha_in,
                          G3MapVectorDouble delta_in,
                          MapCoordReference coord_sys,
                          G3MapVectorDouble & alpha_out,
                          G3MapVectorDouble & delta_out,
                          double vel, double alpha, double delta
                          );

