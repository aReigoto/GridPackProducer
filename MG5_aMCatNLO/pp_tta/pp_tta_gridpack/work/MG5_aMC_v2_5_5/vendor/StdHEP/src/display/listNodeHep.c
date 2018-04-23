/*******************************************************************************
*									       *
* ListNodeHep.c -- Module to display the attributes ( momentum, pt.. ) of a    *
* list of particles characterizing a "nodehep" widget defined in the DispTree  *
* module.  This module must be called from DispTree.			       *
*									       *
* Copyright (c) 1991 Universities Research Association, Inc.		       *
* All rights reserved.							       *
* 									       *
* This material resulted from work developed under a Government Contract and   *
* is subject to the following license:  The Government retains a paid-up,      *
* nonexclusive, irrevocable worldwide license to reproduce, prepare derivative *
* works, perform publicly and display publicly by or for the Government,       *
* including the right to distribute to other Government contractors.  Neither  *
* the United States nor the United States Department of Energy, nor any of     *
* their employees, makes any warrenty, express or implied, or assumes any      *
* legal liability or responsibility for the accuracy, completeness, or         *
* usefulness of any information, apparatus, product, or process disclosed, or  *
* represents that its use would not infringe privately owned rights.           *
*                                        				       *
* Fermilab Nirvana GUI Library						       *
* November 5, 1991							       *
*									       *
* Written by Paul Lebrun & Mark Edel						       *
*									       *
*******************************************************************************/
/*
 * Generated by the ICS builderXcessory (BX).
 *
 *
 * Builder Xcessory 1.0.1.
 *
 */

/*
 * REQUIRED MOTIF INCLUDE FILES
 */
#include <Xm/Xm.h>
#include <X11/Shell.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h> 
#include <Xm/ArrowB.h>
#include <Xm/ArrowBG.h>
#include <Xm/BulletinB.h>
#include <Xm/CascadeB.h>
#include <Xm/CascadeBG.h>
#include <Xm/Command.h>
#include <Xm/CutPaste.h>
#include <Xm/DialogS.h>
#include <Xm/DrawingA.h>
#include <Xm/DrawnB.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/LabelG.h>
#include <Xm/List.h>
#include <Xm/MainW.h>
#include <Xm/MenuShell.h>
#include <Xm/MessageB.h>
#include <Xm/PanedW.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/RowColumn.h>
#include <Xm/Scale.h>
#include <Xm/ScrollBar.h>
#include <Xm/ScrolledW.h>
#include <Xm/SelectioB.h>
#include <Xm/SeparatoG.h>
#include <Xm/Separator.h>
#include <Xm/Text.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <sys/param.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spin/Spin.h"
#include "phase.h"
#include "space.h"
#include "phaseP.h"
#include "dispTree.h"
#include "drawEvent.h"
#include "listNodeHep.h"

/*
** Code
*/
void ListNodeHep(Widget parent, nodehep *ptr)
{
	
    Arg    	args[512];
    int    	argcnt;
    Widget 	retval;
    XmString	xmstr[32];
    XmString	st1;
    Widget	bulletinBoard;
    Widget	form, button;
    Widget	scrolledWindow1;
    Widget	text;
    Widget 	dialogShell;
    char 	*ctext, *cbegin;
    
    int i, j, ip, l, id;
    int *ml;
    double pmom, pt;
    double prm, rm, ppx, ppy, ppz, mass;
    StdHepWindow *w1 = ptr->window;
    PhaseEvent e1;
    PhaseParticle *p1, *ptemp;
    int nlines;
    int lctext;
    char ctmp[80];
    
    
    if (ptr == NULL) {
    	if (w1->nodeWindowShell)
    	    XtDestroyWidget(XtParent(w1->nodeWindowShell));
    	return;
    }
    e1 = w1->event;
    p1 = e1.particles;
    if ((ptr->multiplicity) < 1) return;
    
/* Normal Tree, real event */

    if (w1->modetreedisp == TREEDISPREAL) {
      cbegin = (char *) malloc((ptr->multiplicity)*250 + 100);
      nlines = 5 + 4 * (ptr->multiplicity);
      if (nlines > 30) nlines = 30;
      ctext = cbegin;
      ml = ptr->stdIndex;
      sprintf(ctext, " Kinematics for %d  ", ptr->multiplicity);
      ctext = ctext + strlen(ctext);
      lctext = strlen(ctext);
      hepnam_(&ptr->stdid, ctext, lctext);
      ctext = ctext + 19 ;
      sprintf(ctext,"\n ");
      ctext = ctext + strlen(ctext);
      sprintf(ctext," Standard identity = %d \n ", ptr->stdid);
      ctext = ctext + strlen(ctext);
      ml = ptr->stdIndex;
      ip = *ml;
      ptemp = p1; for(l=0; l< ip; l++, ptemp ++);
       mass =  ptemp->mass;
      sprintf(ctext," Particle(s) Mass = %10.4g \n ", mass);
      ctext = ctext + strlen(ctext);

      ml = ptr->stdIndex;
      for (i=0; i< ptr->multiplicity; i++, ml++){
        sprintf(ctext,"\n");
        ctext = ctext + strlen(ctext);
        ip = *ml;
        ptemp = p1; for(l=0; l< ip; l++, ptemp ++);
      /*
      *place, px, py, pz, 
      * Pt and Rapidity
      */
        ppx = 1.0 * ptemp->px; ppy = 1.* ptemp->py;
        ppz = 1.0 * ptemp->pz;
        pmom = ParticleMomentum(ppx, ppy, ppz);
        pt = ParticlePT(ppx, ppy);
        prm = ParticlePseudorapidity(ppx, ppy, ppz);
        rm = ParticleRapidity(ppx, ppy, ppz, mass);
        sprintf(ctext," Item %d, Index in HEPEVT %4d \n ", i, ip);
        ctext = ctext + strlen(ctext);
      sprintf(ctext," Px, Py, Pz = %10.4g %10.4g %10.4g \n ",
      	ptemp->px, ptemp->py, ptemp->pz);
        ctext = ctext + strlen(ctext);
        sprintf(ctext," Momentum, Pt  = %10.4g %10.4g \n ", pmom, pt);
        ctext = ctext + strlen(ctext);
        sprintf(ctext," Rapidity, Pseudo Rap. = %8.4g %8.4g \n ", 
      				rm, prm);
        ctext = ctext + strlen(ctext);
       }
       if (w1->nodeWindowShell !=NULL)
           XtDestroyWidget(XtParent(w1->nodeWindowShell));
        argcnt = 0;
        XtSetArg(args[argcnt], XmNautoUnmanage, False); argcnt++;
        form = XmCreateFormDialog(parent, "HEP Node Content", args, argcnt);
        w1->nodeWindowShell = form;
      } else {
/*
* Very short documentation about Color code. 
*/
      cbegin = (char *) malloc(500); 
             /* A stupid gues.. Looks verbose enough...*/
      ctext = cbegin;
      ml = ptr->stdIndex;
      sprintf(ctext, " ColorCode Brief Explanation  \n");
      ctext = ctext + strlen(ctext);
      sprintf(ctmp," \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp, " Particle Name is    ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      lctext = strlen(ctmp);
      hepnam_(&ptr->stdid, ctmp, lctext);
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      strcpy(ctext, " \n "); ctext = ctext + strlen(ctext);
      sprintf(ctmp," Standard identity = %d \n ", ptr->stdid);
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      nlines = 10;
      
      id = ptr->stdid; if (id < 0) id = -1 * id;
      
      if ((id == 47) || (id == 25) || (id == 33)) {
      sprintf(ctmp," Supersymetric, Higgs or other \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp,"   undiscovered particles are black \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      
      } else if ((id == 23) || (id == 24) || (id == 9) || (id  == 21)) {
      sprintf(ctmp," Weak Gauge bosons and Gluons are \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp,"  an arbitrary mixture of all 3 RGB colors \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);

      } else if (id == 22)  {
      sprintf(ctmp," The Photon is black. \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The color of a particle is not directly \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," to it's mass or virtuality \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      
      } else if ((id <= 20) && (id > 10)) {
      sprintf(ctmp," Leptons are pure RGB colors ( 3 families only) \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The charge is not represented by a color \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Neutrino are grey. \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      
      } else if (id <= 7)  {
      sprintf(ctmp," The quarks are a mixture of two RGB Colors, \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The color of a particle is not directly \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," to it's mass or virtuality \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The charge is not represented by a color \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      
      } else  {
      sprintf(ctmp," Hadrons:  Color is set by the heaviest quark \n " );
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The charge is not represented by a color \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," Note : The color of a particle is not directly \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      sprintf(ctmp," to it's mass or virtuality \n ");
      strcpy(ctext, ctmp); ctext = ctext + strlen(ctmp);
      }
       if (w1->colorWindowShell !=NULL)
           XtDestroyWidget(XtParent(w1->colorWindowShell));
        argcnt = 0;
        XtSetArg(args[argcnt], XmNautoUnmanage, False); argcnt++;
        form = XmCreateFormDialog(parent, "HEP Node Content", args, argcnt);
        w1->colorWindowShell = form;
     }
      
    if (w1->nodeWindowShell !=NULL)
       XtDestroyWidget(XtParent(w1->nodeWindowShell));
    argcnt = 0;
    XtSetArg(args[argcnt], XmNautoUnmanage, False); argcnt++;
    form = XmCreateFormDialog(parent, "HEP Node Content", args, argcnt);
    w1->nodeWindowShell = form;
    argcnt = 0;
    XtSetArg(args[argcnt], XmNrows, nlines);  argcnt++;
    XtSetArg(args[argcnt], XmNcolumns, 50);  argcnt++;
    XtSetArg(args[argcnt], XmNresizeHeight, False);  argcnt++;
    XtSetArg(args[argcnt], XmNtraversalOn, False); argcnt++;
    XtSetArg(args[argcnt], XmNwordWrap, True);  argcnt++;
    XtSetArg(args[argcnt], XmNscrollHorizontal, False);  argcnt++;
    XtSetArg(args[argcnt], XmNeditMode, XmMULTI_LINE_EDIT);  argcnt++;
    XtSetArg(args[argcnt], XmNeditable, False);  argcnt++;
    XtSetArg(args[argcnt], XmNvalue, cbegin);  argcnt++;
    XtSetArg(args[argcnt], XmNtopAttachment, XmATTACH_FORM);  argcnt++;
    XtSetArg(args[argcnt], XmNleftAttachment, XmATTACH_FORM);  argcnt++;
    XtSetArg(args[argcnt], XmNbottomAttachment, XmATTACH_WIDGET);  argcnt++;
    XtSetArg(args[argcnt], XmNrightAttachment, XmATTACH_FORM);  argcnt++;
    if (nlines < 12) text = XmCreateText(form, "text", args, argcnt);
    else text = XmCreateScrolledText(form, "text", args, argcnt);
    XtManageChild(text);
    
    SET_ONE_RSRC(XtParent(form), XmNtitle, "HEP Node Content");
    XtManageChild(form);
    free (cbegin);
           
}
