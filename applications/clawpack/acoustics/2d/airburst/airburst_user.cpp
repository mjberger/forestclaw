/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "airburst_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#include "../rp/acoustics_user_fort.h"


static
void cb_airburst_output_ascii (fclaw2d_domain_t * domain,
                                fclaw2d_patch_t * this_patch,
                                int this_block_idx, int this_patch_idx,
                               void *user);

void airburst_link_solvers(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *vt = fclaw2d_vt();

    vt->problem_setup = &airburst_problem_setup;  /* Version-independent */
    const user_options_t* user = airburst_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *claw46_vt = fc2d_clawpack46_vt();
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

        claw46_vt->fort_rpn2 = &CLAWPACK46_RPN2;
        claw46_vt->fort_rpt2 = &CLAWPACK46_RPT2;

        claw46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        claw46_vt->fort_bc2       = &CLAWPACK46_BC2;
        claw46_vt->fort_src2      = &CLAWPACK46_SRC2;
        claw46_vt->fort_setaux    = &CLAWPACK46_SETAUX;        

        /* output functions */
        clawpatch_vt->cb_output_ascii  = cb_airburst_output_ascii; 

    }
    else if (user->claw_version == 5)
    {

#if 0
        fc2d_clawpack5_vtable_t    *claw5_vt = fc2d_clawpack5_vt();

        claw5_vt->fort_qinit     = &CLAWPACK5_QINIT;

        claw5_vt->fort_rpn2 = &CLAWPACK5_RPN2;
        claw5_vt->fort_rpt2 = &CLAWPACK5_RPT2;
        }
        else if (user->example == 1)
        {
            fclaw2d_patch_vtable_t  *patch_vt = fclaw2d_patch_vt();
            fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
            
            patch_vt->setup = &airburst_patch_setup;

            claw5_vt->fort_rpn2  = &CLAWPACK5_RPN2_MANIFOLD;
            claw5_vt->fort_rpt2  = &CLAWPACK5_RPT2_MANIFOLD;

            /* Avoid tagging block corners in 5 patch example*/
            clawpatch_vt->fort_tag4refinement = &CLAWPACK5_TAG4REFINEMENT;
            clawpatch_vt->fort_tag4coarsening = &CLAWPACK5_TAG4COARSENING;
        }
#endif
    }

}


void cb_airburst_output_ascii (fclaw2d_domain_t * domain,
                                fclaw2d_patch_t * this_patch,
                                int this_block_idx, int this_patch_idx,
                                void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
    const fclaw_options_t *gparms = fclaw2d_get_options(glob);

    int patch_num;
    int level;
    int mx,my,mbc,meqn;
    double xlower,ylower,dx,dy;
    double *q;
    // char fname[11];

    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);
    
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    char fname[BUFSIZ];
    snprintf (fname, BUFSIZ, "%s.q%04d", gparms->prefix, iframe);
    /* sprintf(fname,"fort.q%04d",iframe); */

    /* The fort routine is defined by a clawpack solver and handles 
       the layout of q in memory (i,j,m) or (m,i,j), etc */
    clawpatch_vt->fort_output_ascii(fname,&mx,&my,&meqn,&mbc,&xlower,&ylower,&dx,&dy,q,
                                    &patch_num,&level,&this_block_idx,
                                    &glob->domain->mpirank);

    
    double yupper = ylower + dx*(mx+1);
    if (yupper > gparms->by)
    {
      double *aux;
      int maux;
      int mitot = mx + 2*mbc;
      int mjtot = my + 2*mbc;
      double time = glob->curr_time;
      double ax = gparms->ax;
      fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
      COLAUX(aux,&maux,&mitot,&mjtot, &xlower, &dx,&time);      
    }


    
}


void airburst_problem_setup(fclaw2d_global_t* glob)
{
    user_options_t* user = airburst_get_options(glob);

    /* rho, bulk are inputs; cc and zz are outputs.  Results are
       stored in a common block */
    AIRBURST_SETPROB(&user->rho, &user->bulk);
}

#if 0
void airburst_patch_setup(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux,*xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    double *xnormals,*ynormals,*xtangents,*ytangents;
    double *surfnormals,*edgelengths,*curvature;

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* Mapped info is needed only for an update */
        return;
    }

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    fclaw2d_clawpatch_metric_data2(glob,this_patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,&edgelengths,
                                   &curvature);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    AIRBURST_SETAUX_MANIFOLD(&mbc,&mx,&my,&xlower,&ylower,
                          &dx,&dy,&maux,aux,
                          xnormals,ynormals,edgelengths,area);
}
#endif
