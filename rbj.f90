MODULE couplage
    IMPLICIT NONE
   rEAL*8,dimension(:,:),allocatable :: Coef_Asperites1d_1,Coef_Asperites1d_2
END MODULE couplage

MODULE ecoulement1d
    IMPLICIT NONE
    real*8,dimension (:),allocatable:: hauteur_old,hauteur,qpsrain1d,HplusZ,widthprime,widthmoyen,Rh,widthmoyen_Maille,dperm,dwidthprime,hauteur_iteration
    real*8,dimension (:),allocatable:: echange1d
    real*8,dimension (:,:),allocatable:: Dechange1d
    integer,dimension (:,:),allocatable:: noeud_inflow,maille_inflow,noeud_outflow,maille_outflow
    integer,dimension (:),allocatable:: Number_inflow,Number_outflow,amont,aval
    real*8,dimension (:,:),allocatable::grad_HplusZ,KP,flux1d
    real*8,dimension (:,:,:),allocatable::dKP,dgradhpz
    real*8 sommepluie1d
END MODULE ecoulement1d

module resolution1d
implicit none
integer iter, itertot
real*8,dimension (:),allocatable:: ax,fd,contrf1d,contrd1d,pdelta,deltaold,calc_deltaold
integer,dimension (:),allocatable:: ai
real*8 omega
real*8,dimension (:),allocatable:: contrf1DTG,contrf1DTD
real*8,dimension (:,:),allocatable::contrhd1d,ContrUP1dTG,ContrDOWN1dTG,ContrUP1dTD,ContrDOWN1dTD
end module

MODULE temporal
    IMPLICIT NONE
    real*8 temps,dt,dtold,start,finish
END MODULE temporal

MODULE maillage_1D
    IMPLICIT NONE
    integer nn1d, nm1d,maxco,nsources,nexits,nb_v1d_m,ntimezonesSOURCES
    real*8,dimension (:,:),allocatable:: cn,cosinusdepente2
    real*8,dimension (:),allocatable:: lengthM_1D,altitude_m1d,lengthM_1DX,distance_exu
    integer,dimension (:,:),allocatable:: mn1d,nconnected,mconnected
    integer,dimension (:),allocatable:: nconnect
    logical interpolerlasource
END MODULE maillage_1D

MODULE Riviere
    IMPLICIT NONE
    real*8,dimension (:,:),allocatable::valimposeZONE
    real*8,dimension (:),allocatable::timeofzonesource,slopesource
    real*8 diflevel,condKP,hds1d,hos1d,hstot1d
    integer,dimension (:),allocatable:: Noeudsource,Noeudsortie,N_NOEUDSinZone
    integer,dimension (:,:),allocatable::NOEUDinZone
    real*8,dimension (:),allocatable:: ValeurIMPOSEE,zonewidth,zonepentel,zonepenter,penter,pentel,width,nman,EpaissSEDIM,KSSEDIM,zonenman,zoneEpaissSEDIM,zoneKSSEDIM
    character*1,dimension (:),allocatable:: GrandeurIMPOSEE   
    integer nwidthzones 
END MODULE Riviere


MODULE Interpolationpluie1d
    IMPLICIT NONE
    integer ntimezonespluie1d, nconditionspluie1d
    logical InterpolerLapluie1d
    REAL*8,dimension(:,:),allocatable :: QpsRAINZONET1d
    integer,dimension(:,:),allocatable :: nmpluvieuses1d
    integer,dimension(:,:,:),allocatable :: maillepluvieuse1d
    rEAL*8,dimension(:),allocatable :: QpsRAINinterp1d,slopeQpsRAIN1d
    rEAL*8,dimension(:),allocatable :: Timeofzonepluie1d
END MODULE Interpolationpluie1d

MODULE Observation1d
    IMPLICIT NONE
   real*8,dimension (:),allocatable:: obst,lengthCANALamont
   integer nobst,nextobs,nbr_noeuds_obs,ncannaux
   integer,dimension(:),allocatable :: noeud_obs,nbr_mailles_cannal,canal,canal_m
   integer,dimension(:,:),allocatable :: numero_maille_cannal
END MODULE Observation1d

MODULE NewtonP
    IMPLICIT NONE
real*8 tmax, epsilonA, epsilonR, Augmentdt, reducdt1, reducdt2,dtmin,dtmax,INIaleatoire
integer mmax,lowlimIT,highlimIT,newtonit
END MODULE NewtonP

MODULE Initialisation1d
    IMPLICIT NONE
    character*15 fichier
    character*12 param_ini1D
    real*8 valueini1D
    integer NexceptionsINI1D
END MODULE Initialisation1d

MODULE bilanmasse1d
    IMPLICIT NONE
     rEAL*8 volume1d_ini,volume1d
     rEAL*8,dimension(:),allocatable :: bilanq1d,bilanQpsRAIN1d,bilanqsourceneumann
END MODULE bilanmasse1d

MODULE Autres
    IMPLICIT NONE
    INTEGER ::  nb_periodes_initial, nb_periodes_final, nb_periodes_max, nb_periodes_sec, nb_periodes            
    REAL :: temps_elapsed
END MODULE Autres


!-----------------------------------------------------------------------------------------------------------------------
!Programme principal
!appel aux différents modules
use temporal
use ecoulement1d
use resolution1d
use riviere
use maillage_1d
use interpolationpluie1d
use observation1d
use newtonp
use bilanmasse1d
use Initialisation1d
implicit none
integer i985,i4155
real*8 variableTEMPO
!Lecture de certains scalaires, donc ceux permettant l'allocation des vecteurs  
call lect0_1D
!allocation de la dimension des vecteurs allocatable
call allocation
!Lecture des vecteurs
!sommepluie1d=0.d0
call lect1_1D
call param_geo1D
call comp_nb_v1d
!initialisation : on calcule la valeur de la source au temps t=0, puis on fait l'initialisation sur tout le domaine
nextobs=1
temps=0.d0
omega=0.d0
bilanq1d=0.d0
bilanqsourceneumann=0.d0
call initialisation_hauteur1d
call calc_volume1d
call plotting
!-----------------------debut du programme
!Boucle en temps
do while (temps<tmax)
                temps=temps+dt 
	                 hauteur_old(1:nn1d)=hauteur(1:nn1D)	      
	                   
	                 !Point de retour en cas de non convergence après mmax iterations             
	                1  continue	                
	                calc_deltaold=0.d0
                    !Initialisation de Picard/Newton
                    NewtonIT=0
                    hauteur(1:nn1D)=hauteur_old(1:nn1D)
                    !calcul du terme puits-source, et de la source
                    call interpolpluie1d
                    call interpolSOURCE
                !Boucle de Newton
                    do                                         
                                    fd=0.d0 ;   Pdelta=0.d0 ;ax=0.d0 ; ai=0.d0 ;contrd1d=0.d0;contrhd1d=0.d0;contrf1d=0.d0    
                                    ContrUP1dTG=0.d0;ContrDOWN1dTG=0.d0;ContrUP1dTD=0.d0;ContrDOWN1dTD=0.d0;contrf1DTG=0.d0;contrf1DTD=0.d0                               
                                    call form_matrice1D         
                                    iter=1 ; itertot=1    

!                                        if (temps.ge.63.75d0) then       
!                                        do i4155=1 ,nb_v1d_m                                    
!                                         write(*,'(a4,i7,a4,i7,a4, f16.10,a4, f16.10)') "riv", ai(i4155), "riv",ai(i4155+nb_v1d_m),"ax",ax(i4155),"fd",fd(ai(i4155))
!                                        enddo
!                                        pause 
!                                        endif
!                       print*, "matrice riviere :"
!do i4155=1 ,nb_v1d_m
!write(*,'(a4,i7,a4,i7, f20.10,a4, f20.10)')"riv", ai(i4155), "riv",ai(i4155+nb_v1d_m),ax(i4155),"fd",fd(ai(i4155))
!if ((isnan(ax(i4155))).or.(isnan(fd(ai(i4155))))) then
!print*, "ohlala noeud", ai(i4155), ax(i4155),fd(ai(i4155))
!pause
!endif
!enddo   

                                    call ndmain(ax,ai,fd,Pdelta,nn1d,nb_v1d_m,iter,itertot) 
                                    do i985=1,nn1d          
                                                hauteur(i985)=pdelta(i985)+hauteur(i985)
                                               if (hauteur(i985).lt.0.d0) then 
                                                   variableTEMPO=hauteur(i985)-pdelta(i985) 
                                                   hauteur(i985)=0.d0
                                                   pdelta(i985)=variableTEMPO                       
                                               endif                                         
                                    end do                                    
                                   
                                    calc_deltaold=calc_deltaold+pdelta
                                     write(*,'(a8,f10.5,a9,i3,a4,f8.5,a8,f6.3,a9,f9.5,a11,f10.5,a7,i6,a3,f12.5,a8,f10.5,a10,f11.5,a10,f11.5)') "temps",temps,"newtonIT", newtonIT," dt ",dt,"omega",omega,"deltamax",MaxVal(dabs(pdelta(1:nn1d)))
                                    newtonIT=newtonIT+1 
                        
                                    if (omega.ne.1.0d0) then
                                        if (MaxVal(dabs(pdelta(1:nn1d))).ge.epsilonA + epsilonr*MinVal(dabs(hauteur(1:nn1d))) ) Then
                                                                                                                    if (newtonIT.ge.mMAX) then 
                                                                                                                                          call GestionDTNONCONV
                                                                                                                                          goto 1
                                                                                                                    endif
                                                                                                                                                          
                                        else 
                                                       exit  
                                        end if
                                    else if (omega.eq.1.0d0) then
                                        exit
                                    endif
                    end do
                !Fin de la boucle de Newton
                omega=0.d0
                Deltaold(1:nn1d)= calc_deltaold(1:nn1d) 
                dtold=dt   
                call calc_volume1d 
                call calc_flux1d    
                call plotting
                call Bilan_masse1d
                call GestionDT
!                
  enddo
!!-----------------------fin du programme
call hauteurs_finales
call temps_simu

end

!--------------------------------------------------------------------------------------------------------
subroutine Gestion_Asperites_1D
use couplage, only : Coef_Asperites1d_1,Coef_Asperites1d_2
use ecoulement1d, only : hauteur_iteration,amont,hauteur_old
use riviere, only : hstot1d,hds1d
use maillage_1d, only : nconnect, nconnected,nn1d,mconnected
implicit none
real*8 hauteurRIV,fraction1,fraction2,fraction3,porosite_asperites
integer i91,j91


do i91=1,nn1d
    do j91=1,nconnect(i91)
        hauteurRIV=hauteur_iteration(amont(mconnected(j91,i91)))  
        !réduction du k écoulement     
        fraction1=  dmin1(1.0d0,(dmax1(hauteurRIV-hds1d,0.d0))/(hstot1d-hds1d))
        Coef_Asperites1d_1(j91,i91)=(fraction1)**(2.0d0*(1.0d0-(fraction1)))        
       if (hstot1d.eq.hds1d)  Coef_Asperites1d_1(j91,i91)=1.0d0
       
       ! réduction du k echange avec le souterrain
        fraction2=  dmin1(1.0d0,(dmax1(hauteurRIV,0.d0))/(hstot1d))
        Coef_Asperites1d_2(j91,i91)=(fraction2)**(2.0d0*(1.0d0-(fraction2)))        
       if ((hstot1d.eq.0.d0).and.(hauteurriv.gt.0.d0) ) then
       Coef_Asperites1d_2(j91,i91)=1.0d0
       elseif ((hstot1d.eq.0.d0).and.(hauteurriv.le.0.d0) ) then
       Coef_Asperites1d_2(j91,i91)=0.0d0
       endif       
    enddo  
enddo
return
end subroutine

!--------------------------------------------------------------------------------------------------------
       subroutine lect0_1D
use temporal
use riviere
use maillage_1d
use interpolationpluie1d
use observation1d
use newtonp
use Initialisation1d
       implicit none
       integer ind

    10 format(a15)
	write(*,*) 'nom du fichier maillage ? '
	read(*,10) fichier
    finish=0.d0
    start=0.d0
    temps=0.d0
call temps_simu
	ind=index(fichier,' ')-1
	
        open(149,file='Input/'//fichier(1:ind)//'/'//"sources1D.txt")
	    open(999,file='Input/'//fichier(1:ind)//'/'//"maillage1d.txt")
	    open(998,file='Input/'//fichier(1:ind)//'/'//"para_1d.txt")
	    open(997,file='Input/'//fichier(1:ind)//'/'//"pluie1d.txt")
	    open(996,file='Input/'//fichier(1:ind)//'/'//"Section1d.txt")
	    open(995,file='Input/'//fichier(1:ind)//'/'//"Tobs.txt")
        open(994,file='Input/'//fichier(1:ind)//'/'//"NewtonParameters.txt")
        open(991,file='Input/'//fichier(1:ind)//'/'//"Cannaux1D.txt")
        open(992,file='Input/'//fichier(1:ind)//'/'//"Noeudobs1D.txt")
        open(993,file='Input/'//fichier(1:ind)//'/'//"Ini1D.txt")
  
 
        open(67,file='Output/'//fichier(1:ind)//'.time')
        open(1950,file='Output/'//fichier(1:ind)//'_flux1d')	
        open(194,file='Output/'//fichier(1:ind)//'_BMasse')
        open(191,file='Output/'//fichier(1:ind)//'_hauteur')   
        open(192,file='Output/'//fichier(1:ind)//'_hauteur.ALL') 

    
    read(991,*) ncannaux
    read(992,*) nbr_noeuds_obs
    print*, "nbr_noeuds_obs",nbr_noeuds_obs
    read(999,*) nn1d,nm1d
    print*, "nn1d",nn1d,"nm1d",nm1d        
    read(998,*) maxco,diflevel
    print*, "maxco",maxco,"diflevel",diflevel
    read (998,*) hds1d
    print*, "hds1d", hds1d
    read (998,*) hos1D
    print*, "hos1d", hos1d
    hstot1d=hds1d+hos1d
    read(998,*) condKP
    print*, "condkp",condkp
    read(998,*) Nexits
    print*,"Nombre de sorties", Nexits
    read(149,*) NtimezonesSources,Nsources,InterpolerlaSource 
    print*, "ntimezonessources", ntimezonessources
    print*, "ns", nsources
    print*, "interpsource", interpolerlasource           
    print*, "nombre de sources", Nsources
    read(997,*) ntimezonespluie1d, nconditionspluie1d,InterpolerLapluie1d
	print*, "ntimezonespluie1d", ntimezonespluie1d
	print*,"nconditionspluie1d",nconditionspluie1d
	print*,"interpolerlapluie1d",interpolerlapluie1d
    read(996,*) nwidthzones
    print*, "nwidthzones",nwidthzones          
	read(995,*)nobst
	print*, "nobst", nobst	
	read(994,*)dt, tmax, epsilonA, epsilonR, Augmentdt, reducdt1, reducdt2,dtmin,dtmax,mmax,lowlimIT,highlimIT,INIaleatoire
	print*, "dt", dt
	print*,"tmax", tmax
	print*,"epsilona", epsilonA
	print*, "epsilonr",epsilonR
	print*, "augmentdt",Augmentdt
	print*,"reducdt1", reducdt1
	print*, "reducdt2",reducdt2
	print*, "dtmin",dtmin
	print*,"dtmax",dtmax
	print*,"mmax",mmax
	print*,"lowlimIT"
	print*,"lowlimIT",lowlimIT
	print*,"highlimIT",highlimIT	
	read (993,*) param_ini1D
	print*, "param_ini", param_ini1D
	read (993,*) valueINi1D
	print*, "valueini", valueini1D
	read (993,*) NexceptionsINI1D
	print*, "nexceptions", nexceptionsini1D

	print*,""
	print*,""
	print*,"*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
	return
	end
!--------------------------------------------------------------------------------------------------------
       subroutine lect1_1D
    use temporal
    use ecoulement1d
    use riviere
    use maillage_1d
    use interpolationpluie1d
    use observation1d
    use newtonp
    use Initialisation1d
       implicit none
       integer iaa,i,j,k,i52,i53,i54,canallocal
       integer indicetemps,indicetemps2
       real*8 tempspluie1d,tempsobs,tempssources
       REAL*8,dimension(:),allocatable ::tempssourcespluie1d,tempsobservation2
       
       !Lecture des coordonnées des noeuds, de la matrice maille noeud, et de la matrice noeuds noeuds      
       do i=1,nn1d
          READ(999,*) iaa,cn(1,i),cn(2,i),cn(3,i)
          cn(3,i)=cn(3,i)-diflevel
 !        write(*,'(i10,3f15.8)') iaa,cn(1,i),cn(2,i),cn(3,i)
       end do      
 !      pause 
       do i=1,nm1d
          READ(999,*) iaa,mn1d(1,i),mn1d(2,i)
 !         print*, iaa,mn1d(1,i),mn1d(2,i)
       end do   
 !      pause    
       do i=1,nn1d
          READ(999,*) iaa,nconnect(i),nconnected(1:nconnect(i),i)
 !          print*, iaa,nconnect(i),nconnected(1:nconnect(i),i)
       end do
 !      pause  
 
  ! Calcul de la longueur de chaque maille, et de l'altitude moyenne de chaque maille
 do i=1,nm1d
 lengthM_1D(i)=((cn(1,mn1d(1,(i)))-(cn(1,mn1d(2,(i)))))**2.0d0  + (cn(2,mn1d(1,(i)))-(cn(2,mn1d(2,(i)))))**2.0d0 + (cn(3,mn1d(1,(i)))-(cn(3,mn1d(2,(i)))))**2.0d0)**0.5d0
 lengthM_1DX(i)=((cn(1,mn1d(1,(i)))-(cn(1,mn1d(2,(i)))))**2.0d0  + (cn(2,mn1d(1,(i)))-(cn(2,mn1d(2,(i)))))**2.0d0)**0.5d0
 altitude_m1d(i)=(cn(3,mn1d(1,(i)))+(cn(3,mn1d(2,(i)))))*0.5d0
 enddo
 
 !Computer, pour chaque noeud, la maille qui le lie au j-ème noeud voisin (matrice noeud maille)
 do i=1,nn1d
   do j=1,nconnect(i)
      do k=1,nm1d
         if (mn1d(1,k).eq.i.and.mn1d(2,k).eq.nconnected(j,i)) then
            mconnected(j,i)=k
            goto 151
         else if (mn1d(2,k).eq.i.and.mn1d(1,k).eq.nconnected(j,i)) then
            mconnected(j,i)=k
            goto 151
         end if   
         !message d'erreur si les noeuds voisins considérés ne partagent en fait pas de maille en commun
         if (k.eq.nm1D) then
             print*, "ERREUR Definition noeuds voisins"
             pause
             stop
         endif    
      end do
      151 continue
   end do
end do    
  
        !Informations sur les sources et sorties
if (nsources.gt.0) then
    do i=1,nsources       
        READ(149,*) Noeudsource(i),GrandeurIMPOSEE(i) 
        print*, Noeudsource(i),GrandeurIMPOSEE(i)      
    enddo            
    do i=1,ntimezonessources
        READ(149,*) timeofzoneSOURCE(i)
        do j=1,nsources
         READ(149,*)  valIMPOSEZone(j,i) 
         print*, valIMPOSEZone(j,i)     
        enddo     
    enddo        
endif
     
       if (nexits.gt.0) then
           do i=1,nexits       
             READ(998,*) Noeudsortie(i)  
           enddo 
       endif 
       
       !Lecture des temps d'observation
       do i= 1, nobst
       read(995,*)obst(i)
       enddo
       
       ! Lecture des noeuds observes
       do i=1,nbr_noeuds_obs
       read(992,*) noeud_obs(i)
       enddo
       
       !Lectures des cannaux d'observation
       if (ncannaux.gt.0) then
           do i=1,ncannaux
                read(991,*) iaa, nbr_mailles_cannal(i)
                do j=1,nbr_mailles_cannal(i)
                read(991,*)    numero_maille_cannal(j,i)
                enddo       
           enddo
       endif 
       !attribution d'un canal à chaque neoud     
       do i52=1,nn1d     
                if (ncannaux.gt.0) then
                    do i53=1,ncannaux
                        do i54=1,nbr_mailles_cannal(i53)
                            if ((mn1d(1,numero_maille_cannal(i54,i53)).eq.i52).or.(mn1d(2,numero_maille_cannal(i54,i53)).eq.i52))then
                                canal(i52)= i53
                                goto 752
                            endif
                        enddo
                    enddo
                else
                canal(i52)= 0                                       
                endif
        752 continue
        enddo
        !attribution d'un canal à chaque maille        
           do i52=1,nm1d     
                if (ncannaux.gt.0) then
                    do i53=1,ncannaux
                        do i54=1,nbr_mailles_cannal(i53)
                            if (numero_maille_cannal(i54,i53).eq.i52)then
                                canal_M(i52)= i53
                                goto 753
                            endif
                        enddo
                    enddo
                else
                canal_m(i52)= 0                                       
                endif
        753 continue
        enddo   

        ! Calcul, pour chaque noeud, de la longueur de canal en amont.
        lengthcanalamont(1:nn1d)=0.d0
        if (ncannaux.gt.0) then
            do i=1,nn1d
                canallocal=canal(i)
                 do j=1,nm1d
                    if  ((altitude_m1d(j).gt.cn(3,i)).and. (canal_m(j).eq.canallocal) ) then
                    lengthCANALamont(i)=lengthcanalamont(i)+lengthM_1D(j)
                    endif
                enddo
            enddo        
        endif
      
       !Lecture des caractéristiques rivière 
       do i= 1, nwidthzones
           read(996,*)zonewidth(i),zonepentel(i),zonepenter(i) , zoneKSSEDIM(i) ,zoneEpaissSEDIM(i),zonenman(i),N_NOEUDSinZone(i)
!           print*,zonewidth(i)
!           print*,zonepentel(i)
!           print*,zonepenter(i)
!           print*, zoneKSSEDIM(i) 
!           print*,zoneEpaissSEDIM(i)
!           print*,zonenman(i)
!           print*,N_NOEUDSinZone(i)
!           pause
            do j=1,N_NOEUDSinZone(i)
               read(996,*)   NOEUDinZone(i,j)
                !print*,  NOEUDinZone(i,j)
           enddo      
       enddo       
       do i=1,nn1D              
           do j=1,nwidthzones
               do k=1,N_NOEUDSinZone(j)
                if (NOEUDinZone(j,k).eq.i) then
                epaisssedim(i)= zoneEpaissSEDIM(j)
                kssedim(i)= zoneKSSEDIM(j)
                nman(i)= zonenman(j)
                width(i)=zonewidth(j)
                penteL(i)=zonepentel(j)
                penteL(i)=penteL(i)*3.14159265358979323846d0/180.0d0                 
                penteR(i)=zonepenter(j)           
                penteR(i)=penteR(i)*  3.14159265358979323846d0/180.0d0   
                goto 152  
                endif       
               enddo               
               !message d'erreur si le noeud considéré n'est dans aucune zone
               if (j.eq.nwidthzones) then
               print*, "ERREUR Definition zones Rivière"
               pause
               stop
               endif               
           enddo
           152 continue
           !print*, cos(pentel(i)), cos(penter(i)),width(i)
       enddo          
 
       !Lecture de la pluie1d
          if (nconditionspluie1d>0) then
       do i=1,ntimezonespluie1d
            read(997,*) timeofzonepluie1d(i)
 !           print*, timeofzonepluie1d(i)
            do j= 1,nconditionspluie1d
               read(997,*) QpsRAINZONET1d(j,i),nmpluvieuses1d(j,i)
 !               print*, QpsRAINZONET1d(j,i),nmpluvieuses1d(j,i)
               do k=1,nmpluvieuses1d(j,i)
               read(997,*) maillepluvieuse1d(j,i,k)
  !             sommepluie1d=sommepluie1d+lengthM_1Dx(maillepluvieuse1d(j,i,k))
 !              print*, maillepluvieuse1d(j,i,k)
               enddo
  !             print*, sommepluie1d
 !              sommepluie1d=0.d0
    
            enddo
       enddo    
    endif   
    
    
!Nous allons inclure aux temps d'observation 1d les temps correspondants aux changements de pluie1d et de condition limite.
        allocate(tempssourcespluie1d(ntimezonespluie1d+ntimezonessources+nobst))
        !d'abord, on s'occupe de tempspluie1d
        indicetemps=0
        do i=1,ntimezonespluie1d
            indicetemps=indicetemps+1
            tempspluie1d=timeofzonepluie1d(i)
            tempssourcespluie1d(indicetemps)=tempspluie1d
        enddo
        
        !maintenant on ajoute les tempssources différents des tempspluie1d déja computés        
        do i=1,ntimezonessources
            tempssources=timeofzonesource(i)
            do j=1,ntimezonespluie1d
            tempspluie1d=timeofzonepluie1d(j)
            if (tempssources.eq.tempspluie1d )then
                goto 2590
            endif            
            enddo
            indicetemps=indicetemps+1
            tempssourcespluie1d(indicetemps)=tempssources
            2590 continue    
        enddo        
       
        !maintenant on ajoute les temps réels d'observation qui ne sont pas égaux aux temps pluie1d/BC !
        indicetemps2=indicetemps
        do i=1,nobst
            tempsobs=obst(i)
            do j=1,indicetemps2
                if (tempsobs.eq.tempssourcespluie1d(j)) then
                goto 2600
                endif
            enddo       
            indicetemps=indicetemps+1
            tempssourcespluie1d(indicetemps)=tempsobs 
            2600 continue
        enddo
        
        !on met ce vecteur dans un vecteur de la bonne taille
        allocate(Tempsobservation2(indicetemps))
        do i=1,indicetemps
            tempsobservation2(i)=tempssourcespluie1d(i)
        enddo
        deallocate(tempssourcespluie1d)
                
        !Maintenant on classe tempsobservation2 par ordre croissant 
        call sort_REAL8array(tempsobservation2,size(tempsobservation2,1))     
        !puis, on remplace obst !
        deallocate(obst)
        nobst=indicetemps
        allocate(obst(nobst))
        do i=1,nobst
         obst(i)=tempsobservation2(i)
        enddo     
        deallocate(tempsobservation2)
!print*,"pi"
!do i=1,nobst
!print*, obst(i)
!enddo
!print*,"pi"
!pause
	print*,""
	print*,"*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
	
	return
	end


!-----------------------------------------------------------------------------------------------------
 subroutine initialisation_hauteur1d
 use Initialisation1d
 use maillage_1D, only : nn1D,nsources,distance_exu,cn
 use ecoulement1d, only: hauteur
 use riviere, only : grandeurimposee,noeudsource,valeurimposee,noeudsortie
 use newtonp, only: INIaleatoire
implicit none
integer Noeud_exception,i
real*8 valueexception1D,maxD,aleatoire
integer i34,j34,i35

!SI LA HAUTEUR EST IMPOSEE
if ((param_ini1D).eq."Hauteur_lame") then
   hauteur(1:nn1D)=valueini1D
    if (nexceptionsini1D>0) then
        do i34=1,nexceptionsini1D
            read (993,*) Noeud_exception,valueexception1D
            hauteur(Noeud_exception)=valueexception1D
        enddo
    endif
else if ((param_ini1D).eq."Coef_Pente_H") then
    do i34=1,nn1d
        Distance_exu(i34)= ((cn(1,noeudsortie(1))-cn(1,i34))**2.0d0+(cn(2,noeudsortie(1))-cn(2,i34))**2.0d0+((cn(3,noeudsortie(1))-cn(3,i34))**2.0d0))**0.5d0
    enddo
    maxD= maxval(Distance_exu(1:nn1d))
    do i34=1,nn1d
        hauteur(i34)=10.0d-10+Distance_exu(i34)*valueini1d/maxD
    enddo
    if (nexceptionsini1D>0) then
        do i34=1,nexceptionsini1D
            read (993,*) Noeud_exception,valueexception1D
            hauteur(Noeud_exception)=valueexception1D
        enddo
    endif
else
    print*, "ERROR INITIALISATION : check input"
    pause
    stop
endif  

    do i34=1,nn1d
aleatoire= dabs((dcos(cn(1,i34)**2.0d0+cn(2,i34))))*INIaleatoire
hauteur(i34)=hauteur(i34)+aleatoire
enddo


!do i34=1,nn1d
!print*, hauteur(i34)
!enddo
!print*, "ini"
return
 end subroutine
!-----------------------------------------------------------------------------------------------------
 subroutine interpolpluie1d
 use Interpolationpluie1d
 use maillage_1d, only : nm1D,lengthM_1D,lengthM_1DX
 use ecoulement1d, only : qpsrain1d
 use temporal, only : temps
 implicit none
 integer i31,j31,k31,i97
 
qpsrain1d=0.d0
 
 if (nconditionspluie1d>0) then 
        if (InterpolerLapluie1d.eq..TRUE.) then
         do i31=1,ntimezonespluie1d     
             if (i31==ntimezonespluie1d) then
                     do j31=1,nconditionspluie1d
                        QpsRAINinterp1d(j31)=QpsRAINZONET1d(j31,i31)
                        do k31 =1,nmpluvieuses1d(j31,i31)
                           qpsrain1d(maillepluvieuse1d(j31,i31,k31))= QpsRAINinterp1d(j31)
                        enddo 
                      enddo
                      exit
             elseif ((temps >= timeofzonepluie1d(i31)) .and.(temps < timeofzonepluie1d(i31+1)))  then
                 do j31=1,nconditionspluie1d
                     slopeQpsRAIN1d(j31)=(QpsRAINZONET1d(j31,i31+1)-QpsRAINZONET1d(j31,i31))/(timeofzonepluie1d(i31+1)-timeofzonepluie1d(i31))
                     QpsRAINinterp1d(j31)=QpsRAINZONET1d(j31,i31)+slopeQpsRAIN1d(j31)*(temps-timeofzonepluie1d(i31))
                     do k31 =1,nmpluvieuses1d(j31,i31)
                      qpsrain1d(maillepluvieuse1d(j31,i31,k31))= QpsRAINinterp1d(j31)
                     enddo                                      
                 enddo
                 exit
             endif
         enddo 
        else
          do i31=1,ntimezonespluie1d     
             if (i31==ntimezonespluie1d) then
                     do j31=1,nconditionspluie1d
                        QpsRAINinterp1d(j31)=QpsRAINZONET1d(j31,i31)
                        do k31 =1,nmpluvieuses1d(j31,i31)
                           qpsrain1d(maillepluvieuse1d(j31,i31,k31))= QpsRAINinterp1d(j31)
                        enddo 
                      enddo
                      exit
             elseif ((temps >= timeofzonepluie1d(i31)) .and.(temps.le.timeofzonepluie1d(i31+1)))  then
                     do j31=1,nconditionspluie1d
                        QpsRAINinterp1d(j31)=QpsRAINZONET1d(j31,i31)
                        do k31 =1,nmpluvieuses1d(j31,i31)
                           qpsrain1d(maillepluvieuse1d(j31,i31,k31))= QpsRAINinterp1d(j31)
                        enddo 
                      enddo
                      exit
             endif
         enddo
        endif
endif


do j31=1,nm1d
!print*, "maille",j31,"QpsRAINRAW",qpsrain1d(j31)
qpsrain1d(j31)=qpsrain1d(j31)*(lengthM_1DX(j31))/(lengthM_1D(j31))
!print*, "qpsrain1d corrige",qpsrain1d(j31)*lengthM_1D(j31)
enddo

!print*, sum(qpsrain1d(:)*lengthM_1D(:))*60.0d0
!pause
end subroutine
!-----------------------------------------------------------------------------------------------------
       subroutine allocation
       use ecoulement1d
        use temporal
        use resolution1d
        use riviere
        use maillage_1d
        use interpolationpluie1d
        use observation1d
        use newtonp
        use Initialisation1d
        use couplage
        use bilanmasse1d
        implicit none  
    
  allocate (Coef_Asperites1d_1(maxco,nn1d),Coef_Asperites1d_2(maxco,nn1d),dperm(nn1d))
 allocate (dgradhpz(maxco,nn1d,2)) 
 allocate (ContrUP1dTG(maxco,nn1d),ContrDOWN1dTG(maxco,nn1d),ContrUP1dTD(maxco,nn1d),ContrDOWN1dTD(maxco,nn1d),contrf1DTG(nn1d),contrf1DTD(nn1d))       
 allocate (zoneNman(nwidthzones),zoneepaisssedim (nwidthzones),zonekssedim (nwidthzones),nman(nn1d),kssedim(nn1d),epaisssedim(nn1d))        
 allocate (cosinusdepente2(maxco,nn1d))       
 allocate (distance_exu(nn1d))       
 allocate(timeofzonepluie1d(ntimezonespluie1d),maillepluvieuse1d(nconditionspluie1d,ntimezonespluie1d,nm1d),nmpluvieuses1d(nconditionspluie1d,ntimezonespluie1d),QpsRAINZONET1d(nconditionspluie1d,ntimezonespluie1d),QpsRAINinterp1d(nconditionspluie1d),slopeQpsRAIN1d(nconditionspluie1d))
!obst
 allocate (obst(nobst),cn(3,nn1d),mn1d(2,nm1d)) 
 allocate (nconnected(maxco,nn1d),mconnected(maxco,nn1d),nconnect(nn1d))
 allocate ( Noeudsortie(nexits),Noeudsource(nsources),ValeurIMPOSEE(nsources),GrandeurIMPOSEE(nsources)   )
 allocate( zonewidth(nwidthzones),zonepentel(nwidthzones),zonepenter(nwidthzones),  N_NOEUDSinZone(nwidthzones)   )        
 allocate (NOEUDinZone (nwidthzones,nn1D) )
 allocate(hauteur(nn1d),lengthM_1D(nm1d),altitude_m1d(nm1d))
 allocate(qpsrain1d(nm1d),dwidthprime(nn1d))
 allocate (hplusZ(nn1d))
 allocate ( noeud_inflow(maxco,nn1d),maille_inflow(maxco,nn1d),noeud_outflow(maxco,nn1d),maille_outflow(maxco,nn1d),Number_inflow(nn1d),Number_outflow(nn1d)  )
 allocate (widthprime(nn1d),widthmoyen(nn1d),width(nn1d),penter(nn1d),pentel(nn1d) )    
 allocate (Rh(nn1d),widthmoyen_Maille(nm1d),amont(nm1d),aval(nm1d))
 allocate (kp(maxco,nn1d),grad_HplusZ(maxco,nn1d),dkp(maxco,nn1d,2) )
 allocate (hauteur_old(nn1d),hauteur_iteration(nn1d),lengthM_1Dx(nm1d))
 !toute la ligne suivante
 allocate (ax(10*nn1d),ai(20*nn1d),fd(nn1d),contrf1d(nn1d),contrd1d(nn1d),contrhd1d(maxco,nn1d))
 !deltaold,calc_deltaold
 allocate( noeud_obs(nbr_noeuds_obs),deltaold(nn1d),calc_deltaold(nn1d) )
 allocate (pdelta(nn1d),flux1d(maxco,nn1d),bilanq1d(nn1d),bilanQpsRAIN1d(nm1d))
 allocate( nbr_mailles_cannal(ncannaux),numero_maille_cannal(nm1d,ncannaux),canal(nn1d) )
 allocate (lengthCANALamont(nn1d),canal_m(nm1d))
 allocate (valimposezone(nsources,ntimezonessources)   )
 allocate(timeofzonesource(ntimezonessources),slopesource(nsources),bilanqsourceneumann(nn1d))
 ! La taille nn1d est ici provisoire
 allocate (echange1d(nn1d),Dechange1d(2,nn1d))

	return
	end
	
!-----------------------------------------------------------------------------------------------------
 subroutine form_matrice1D  
 use ecoulement1d, only : hauteur, hauteur_iteration,hauteur_old
 use maillage_1d, only : nn1d
 use resolution1d, only : omega
 implicit none
 hauteur_iteration(1:nn1d)=hauteur(1:nn1d)*(1.0d0-omega)+hauteur_old(1:nn1d)*omega
  call param_hydro1d
  call comp_Coefs       
	end

!-----------------------------------------------------------------------------------------------------
subroutine param_hydro1d
 use maillage_1D,only : nconnect,nconnected,mconnected,lengthM_1D,mn1D,nm1D,nn1d,cn
 use ecoulement1d, only : dgradhpz,widthmoyen_Maille,rh,dperm,dwidthprime,hauteur_iteration,amont, aval,kp,dkp,widthPRIME,widthmoyen,HplusZ,grad_HplusZ
 use riviere, only : hds1d,width,nman,pentel,penter,condkp
 use temporal, only : temps
 use newtonp, only : newtonit
 use resolution1d, only: omega
 use couplage, only :  Coef_Asperites1d_1
 implicit none
 integer i,j,noeudamont,noeudaval
 real*8 penteldown,pentelup,penterup,penterdown,zup,zdown,ldown,lup,lgthmail,nmanlocal
 real*8 t2,t4,t6,t10,t12,t14,t15,t17,t18,t20,t21,t22,t24,t26,t28,t29,t31,t32,t33,t35,t36,t37,t38,t50,t61,t64
 real*8 casp,lengthloc,dgradhpzup,dgradhpzdown,hsaufhdsup,hsaufhdsdown,dhsaufhdsup,dhsaufhdsdown,grad,perm
 integer mailleconnec
 
 do i=1,nn1d
      !H plus Z
     HplusZ(i)=cn(3,i)+hauteur_iteration(i) 
       !Lprime
     widthprime(i)=width(i)  + hauteur_iteration(i)* dwidthprime(i)
    ! print*, "widthp",widthprime(i)
     !Lmoyen
     widthmoyen(i)= 0.5d0*(  width(i)+widthprime(i))
      !Rayon Hydraulique
      perm=  width(i) + hauteur_iteration(i)*dperm(i)
     Rh(i)=(hauteur_iteration(i)*widthmoyen(i) )/perm 
    ! print*, "rh",rh(i)
     
enddo
  !Moyenne sur les deux noeuds d'une maille de Lmoyen 
 ! print*, nm1d
 do i=1,nm1d
 !print*, i
    widthmoyen_Maille(i)= 0.5d0*(  widthmoyen(mn1D(1,i))+widthmoyen(mn1D(2,i))  )
 enddo  
 !Détermination du noeud amont pour chaque maille
 do i=1,nm1d
     if (HplusZ(mn1D(1,i)).gt.HplusZ(mn1D(2,i))) then
         amont(i)=mn1D(1,i)
         aval(i)=mn1D(2,i)
     else
         amont(i)=mn1D(2,i) 
         aval(i)=mn1D(1,i)
     endif
 enddo
 !Gradient de H+Z
 do i=1,nn1d
     do j=1,nconnect(i)
     mailleconnec=mconnected(j,i)
     lengthloc=lengthM_1D(mailleconnec)
        grad_HplusZ(j,i)=  -1.0d0*dabs((-HplusZ(i)+HplusZ(nconnected(j,i)))/lengthloc)
        dgradhpz(j,i,1)= -1.0d0/lengthloc     ! Dérivée par rapport au noeud amont
        dgradhpz(j,i,2)= 1.0d0/lengthloc      ! Dérivée par rapport au noeud aval
     enddo
 enddo

call Gestion_Asperites_1D

 !KP et Dérivée de KP par rapport à H : Voir MAPLE
do i=1,nn1d
    do j=1,nconnect(i)
        noeudamont=amont(mconnected(j,i))
        noeudaval=aval(mconnected(j,i))
        penteLDOWN=penteL(noeudaval)
        penteLUP=penteL(noeudamont)
        penterDOWN=penter(noeudaval)
        penterUP=penter(noeudamont)
        LUP=width(noeudamont)
        LDOWN=width(noeudaval)
        ZUP=cn(3,noeudamont)
        ZDOWN=cn(3,noeudaval)
        lgthmail=lengthM_1D(mconnected(j,i)) 
        Nmanlocal=(Nman(noeudamont)+Nman(noeudaval))*0.5d0
        grad= grad_HplusZ(j,i)
        dgradhpzup= dgradhpz(j,i,1)
        dgradhpzdown= dgradhpz(j,i,2)
        hsaufhdsUP=dmax1(hauteur_iteration(noeudamont)-hds1d,0.d0)
        if (hsaufhdsup.eq.0.d0) then
            kp(j,i)=0.d0
            dkp(j,i,1)=0.d0
            dkp(j,i,2)=0.d0
            goto 245
        endif
        
        if (hsaufhdsup.gt.0.d0) then
            dhsaufhdsup =1.0d0
        else
             dhsaufhdsup =0.0d0
        endif
        hsaufhdsDOWN=dmax1(hauteur_iteration(noeudaval)-hds1d,0.d0)
        if (hsaufhdsDOWN.gt.0.d0) then
             dhsaufhdsDOWN=1.0d0
        else
             dhsaufhdsDOWN=0.d0
        endif
      t2 = dtan(penteLDOWN)
      t4 = dtan(penteRDOWN)
      t6 = 0.1D1 / t2 + 0.1D1 / t4
      t10 = dtan(penteLUP)
      t12 = dtan(penteRUP)
      t14 = 0.1D1 / t10 + 0.1D1 / t12
      t15 = hsaufhdsup * t14
      t17 = LDOWN / 0.2D1 + hsaufhdsdown * t6 / 0.4D1 + LUP / 0.2D1 + t15 / 0.4D1
      t18 = t17 * hsaufhdsup
      t20 = LUP + t15 / 0.2D1
      t21 = t20 * hsaufhdsup
      t22 = dsin(penteLUP)
      t24 = dsin(penteRUP)
      t26 = 0.1D1 / t22 + 0.1D1 / t24
      t28 = hsaufhdsup * t26 + LUP
      t29 = 0.1D1 / t28
      t31 = (t21 * t29) ** (0.1D1 / 0.3D1)
      t32 = t31 ** 2
      t33 = 0.1D1 / Nmanlocal
      t35 = dsqrt(-grad)
      t36 = 0.1D1 / t35
      t37 = t32 * t33 * t36
      kp(j,i)= t18 * t37
      t38 = dhsaufhdsup * t14
      t50 = t28 ** 2
      t61 = t18 * t32
      t64 = -0.1D1 / t35 / grad
     
      if (noeudamont.eq.i) then
          dkp(j,i,1) = t38 * hsaufhdsup * t37 / 0.4D1 + t17 * dhsaufhdsup * t37 + 0.2D1 / 0.3D1 * t18 * (t38 * hsaufhdsup * t29 / 0.2D1 + t20 * dhsaufhdsup * t29 - t21 * dhsaufhdsup * t26 / t50) / t31 * t33 * t36 + t61 * dgradhpzup * t33 * t64 / 0.2D1
          dkp(j,i,2) = dhsaufhdsdown * t6 * hsaufhdsup * t37 / 0.4D1 + t61 * dgradhpzdown * t33 * t64 / 0.2D1
    elseif (noeudamont.eq.nconnected(j,i)) then
           dkp(j,i,2) = t38 * hsaufhdsup * t37 / 0.4D1 + t17 * dhsaufhdsup * t37 + 0.2D1 / 0.3D1 * t18 * (t38 * hsaufhdsup * t29 / 0.2D1 + t20 * dhsaufhdsup * t29 - t21 * dhsaufhdsup * t26 / t50) / t31 * t33 * t36 + t61 * dgradhpzup * t33 * t64 / 0.2D1
           dkp(j,i,1) = dhsaufhdsdown * t6 * hsaufhdsup * t37 / 0.4D1 + t61 * dgradhpzdown * t33 * t64 / 0.2D1
       else
                print*, "error in kp or dkp"
                pause
                stop        
       endif 

!print*, "haut_iter",hauteur_iteration(noeudamont),"hds",hds1d
!print*, "i",i,"amont", noeudamont
!print*, t2,t2
!print*, t4,t4
!print*, "t6",t6
!print*, "t10",t10
!print*, "t12",t12
!print*, "t14",t14
!print*, "hsaufhdsup",hsaufhdsup
!print*, "hds",hds1d
!print*, "t15",t15
!print*, "t17",t17
!print*, "t18",t18
!print*, "t20",t20
!print*, "t21",t21
!print*, "t22",t22
!print*, "t24",t24
!print*, "t26",t26
!print*, "t28",t28
!print*, "t29",t29
!print*, "t31",t31
!print*, "t32",t32
!print*, "t33",t33
!print*, "t35",t35
!print*, "t36",t36
!print*, "t37",t37
!print*, "t50",t50
!print*, "t61",t61
!print*, "t64",t64
!print*, "kp",kp(j,i)
!print*, "dkp",dkp(j,i,1),dkp(j,i,2)

if (kp(j,i).gt.condKP) then
kp(j,i)=condKP
dkp(j,i,1)=condKP
dkp(j,i,2)=condKP
endif

!réduction a cause de l'obstruction
casp= Coef_Asperites1d_1(j,i)
kp(j,i)=kp(j,i) * casp
dkp(j,i,1)=dkp(j,i,1) * casp
dkp(j,i,2)=dkp(j,i,2)* casp
245 continue
enddo
enddo 
	end	
!-----------------------------------------------------------------------------------------------------
subroutine param_geo1d
 use maillage_1D,only : nn1d
 use ecoulement1d, only : dwidthprime,dperm
 use riviere, only : pentel,penter
 implicit none
 integer i

 do i=1,nn1d    
     dwidthprime(i)=1.0d0/dtan(pentel(i))+1.0d0/dtan(penteR(i)) 
     dperm(i)=( 1.0d0/dSIN(pentel(i))+1.0d0/dSIN(penteR(i)))  
enddo
	end	
!-----------------------------------------------------------------------------------------------------
subroutine comp_Coefs
use ecoulement1d, only: qpsrain1d,amont, aval,kp,hauteur_old,widthprime,dkp,hauteur_iteration,dwidthprime
use maillage_1d, only : nb_v1d_m,nn1d, mconnected, nconnected,lengthM_1D,cn,nsources,nexits,nconnect
use resolution1d, only : fd, ax, ai,contrd1d,contrhd1d,contrf1d,omega,contrf1dTG,contrf1dTD,ContrUP1dTD,ContrUP1dTG,ContrDOWN1dTD,ContrDOWN1dTG
use riviere, only : hds1d,grandeurimposee, valeurimposee,noeudsortie,noeudsource
use temporal, only : dt,temps
use newtonp, only: newtonit
implicit none
integer i,j,k,cal_ne,noeudamont,noeudaval,i6,cal_nediag
real*8 hauteurdown,hauteurup,zdown,zup,lup,ldown,lgthmail,kplocal,QpsRAINlocal,hauteuroldup,hauteurolddown
real*8 lprimeup,lprimedown,dkplocaldown,dkplocalup,dlprime
real*8 hauteurdiag,hauteurolddiag,lprimediag,dhauteurdiagdown,dhauteurdiagUP,dlprimedown,dlprimeup,fonctionsigne
real*8 t1,t3,t4,t7,t9,t16,t15,t13,t20,t21,t22,t10,t14,t23
 cal_ne=0   
  
 !Calcul, pour chaque couple (i,j), du coefficient à rentrer dans la matrice.
do i=1,nn1d
!TRAITEMENT DES SOURCES A HAUTEUR IMPOSEE
if (nsources.gt.0) then
do k=1,nsources
    if (i.eq.noeudsource(k)) then
        if (grandeurimposee(k).eq."H") then
            contrd1d(i)=1.0d0
            contrf1d(i)=valeurimposee(k)-hauteur_iteration(i)   
            !METTRE LA MATRICE EN FORME CREUSE
            cal_ne = cal_ne+1
            ai(cal_ne) = i
            ai(nb_v1d_m+cal_ne) = i
            ax(cal_ne) = contrd1d(i)
            fd(i) = -contrf1d(i)   
            goto 553   
        endif
    endif    
enddo
endif


!TRAIEMENT DU CAS GENERAL
    do j=1,nconnect(i)    
            noeudamont=amont(mconnected(j,i))
            noeudaval=aval(mconnected(j,i))
            hauteurUP=hauteur_iteration(noeudamont) 
            hauteurDOWN=hauteur_iteration(noeudaval)
            ZUP=cn(3,noeudamont)
            ZDOWN=cn(3,noeudaval)
            lgthmail=lengthM_1D(mconnected(j,i))
            kplocal=kp(j,i)
            QpsRAINlocal=qpsrain1d(mconnected(j,i))
            hauteurdiag=hauteur_iteration(i)
            if (hauteurdiag.eq.0.d0) then
                hauteurolddiag=0.d0
            else
                hauteurolddiag=hauteur_old(i)
            endif
            Lprimediag=widthprime(i)
            do k=1,nsources
                if ((i.eq.noeudsource(k)).and.(grandeurimposee(k).eq."F")) then
!                print*, "i",i,"noeudsource(k)",noeudsource(k)
                QpsRAINlocal=QpsRAINlocal+ 2.0d0*valeurimposee(k)/( Lprimediag*lgthmail ) 
!                print*, "j", j
!                print*, "k",k
!                print*, "qpsrainlocal",qpsrainlocal
!                print*, "valimp", valeurimposee(k)
!                pause
                exit 
                endif
            enddo            

            if (noeudamont.eq.i) then
                fonctionsigne=-1.0d0
                dkplocalup =dkp(j,i,1)
                dkplocaldown=dkp(j,i,2) 
                dhauteurdiagdown=0.d0
                dlprimedown=0.d0
                dhauteurdiagup= 1.0d0
                dlprimeup=dwidthprime(i)
            else
                fonctionsigne=1.0d0
                dkplocalup =dkp(j,i,2)
                dkplocaldown=dkp(j,i,1) 
                dhauteurdiagdown=1.0d0
                dlprimedown=dwidthprime(i)
                dlprimeup=0.d0
                dhauteurdiagup=0.d0
            endif            
      t1 = Lprimediag * lgthmail
      t3 = 0.1D1 / dt
      t4 = (hauteurdiag - hauteurOLDdiag) * t3
      t7 = Fonctionsigne * kplocal
      t9 = 0.1D1 / lgthmail
      t10 = (HauteurDOWN + ZDOWN - HauteurUP - ZUP) * t9
      contrf1DTG(i) =contrf1DTG(i)+ t1 * t4 / 0.2D1 + t7 * t10
      contrf1DTD(i) =contrf1DTD(i) -t1 * QpsRainlocal / 0.2D1
      t14 = dLprimeDOWN * lgthmail
      t22 = t7 * t9
      ContrDOWN1dTG(j,i) = t14 * t4 / 0.2D1 + t1 * dHauteurDiagDOWN * t3 / 0.2D1 + Fonctionsigne * dKplocalDOWN * t10 + t22
      t23 = dLprimeUP * lgthmail
      ContrUP1dTG(j,i) = t23 * t4 / 0.2D1 + t1 * dHauteurDiagUP * t3 / 0.2D1 + Fonctionsigne * dKplocalUP * t10 - t22
      ContrDOWN1dTD(j,i) = -t14 * QpsRainlocal / 0.2D1
      ContrUP1dTD(j,i) = -t23 * QpsRainlocal / 0.2D1   
      call traite_sorties2(i,j,noeudamont,noeudaval,zup,zdown,lgthmail,lprimediag)
      !call traite_sorties2(i,j,noeudamont,noeudaval,hauteurdiag,dhauteurdiagdown,dhauteurdiagup,zup,zdown,lgthmail,lprimediag)
     
!      if ((i.eq.2).and.(noeudamont.ne.i)) then
!      print*, "i", i, "j", nconnected(j,i)
!      print*, "amont", noeudamont,"ch", hauteurup+zup,"haut",hauteurup
!      print*, "aval", noeudaval,"ch", hauteurdown +zdown,"haut",hauteurdown
!      print*,   "ContrUP1dTG", ContrUP1dTG(j,i)
!      print*,  "ContrDOWN1dTG",ContrDOWN1dTG(j,i) 
!      print*,   "ContrUP1dTD", ContrUP1dTD(j,i)
!      print*, "ContrDOWN1dTD", ContrDOWN1dTD(j,i) 
!      print*, "contrf1DTG",contrf1DTG(i)
!      print*, "contrf1DTD",contrf1DTD(i)
!        print*, "t1",t1
!        print*, "t3",t3
!        print*, "t4",t4
!        print*, "t7",t7
!        print*, "t9",t9
!        print*, "t13",t13
!        print*, "t21",t21
!         print*, "t22",t22
!         print*, dhauteurdiagup,"dhauteurdiagup"
!         print*, "test", - dKplocalUP * t7 * t9 + t21
!      endif
        enddo  
        
        !METTRE LA MATRICE EN FORME CREUSE
        cal_ne=cal_ne+1
        cal_neDiag=cal_ne
        ai(cal_neDiag) = i
        ai(nb_v1d_m+cal_neDiag)=i
        fd(i)= -1.0d0* (contrf1DTG(i)+ contrf1DTD(i))
        do j=1,nconnect(i) 
            cal_ne = cal_ne+1
            ai(cal_ne) = i
            ai(nb_v1d_m+cal_ne) = nconnected(j,i)
            if (amont(mconnected(j,i)).eq.i) then
                ax(cal_nediag)=ax(cal_nediag)+ ContrUP1dTG(j,i) + ContrUP1dTD(j,i)
                ax(cal_ne)=ContrDOWN1dTD(j,i) + ContrDOWN1dTG(j,i)
            else
                ax(cal_nediag)=ax(cal_nediag)+ ContrDOWN1dTG(j,i) + ContrDOWN1dTD(j,i)
                ax(cal_ne)=ContrUP1dTD(j,i) + ContrUP1dTG(j,i)        
            endif
        enddo       
           553 continue
 enddo
 end	
!---------------------------------------------------------------------------------------------------
subroutine traite_sorties2(i,j,noeudamont,noeudaval,zup,zdown,lgthmail)
use maillage_1d, only : nexits,mconnected
use ecoulement1d, only : hauteur_iteration,amont,dwidthprime,dperm
use riviere, only :  noeudsortie,nman,width,pentel,penter,hds1d
use resolution1d, only : contrf1dTD,ContrUP1dTD,ContrDOWN1dTD
use couplage, only :  Coef_Asperites1d_1
use temporal, only : temps
 implicit none
 integer i,k,noeudamont,noeudaval,j
 real*8 u1,u2,u3,u5,u6,u7,u8,u9,u10,u12,u22
 real*8 pentelocal,zup,zdown,lgthmail,nmanlocal,dlmoyenup,dlmoyendown,lmoyendiag,perimMouilleDIAG,dpermShds
 real*8 dlp,dpm,hauteurman,dhauteurmandown,dhauteurmanup,permshds,dpermshdsup,dpermshdsdown,casp
       !Lprime
     
do k=1,nexits
    if (i.eq.noeudsortie(k)) then
     casp=Coef_Asperites1d_1(1,i)
     Hauteurman=  (dmax1(hauteur_iteration(i)-hds1d,0.d0))
     Nmanlocal=(Nman(i)) 
     pentelocal=dabs(zup-zdown)/lgthmail
     dlp = dwidthprime(i)
     dpm= dperm(i)
     lmoyendiag= ((width(i)  + hauteurman* dlp )+width(i))/2.0d0
     permshds=  width(i) + dmax1(hauteur_iteration(i)-hds1d,0.d0)*dpm
     if (hauteurman.eq.0.d0) then           
         dHauteurManDOWN=0.d0
         dHauteurManUP=0.d0
     else 
         if (amont(mconnected(j,i)).eq.i ) then
             dHauteurManDOWN=0.d0
             dHauteurManUP=1.d0 
         else
             dHauteurManDOWN=1.d0
             dHauteurManUP=0.d0 
         endif   
     endif
     dpermShdsdown=dhauteurmandown*dpm  
     dpermshdsup=dhauteurmanup*dpm    
     dLmoyenUP= dhauteurmandown*dlp/2.0d0
     dLmoyenDOWN=dhauteurmanup*dlp/2.0d0 
     
      u1 = HauteurMan * Lmoyendiag
      u2 = u1 ** (0.1D1 / 0.3D1)
      u3 = u2 ** 2.0d0
      u5 = dsqrt(pentelocal)
      u6 = u3 * u1 * u5
      u7 = 0.1D1 / Nmanlocal
      u8 = permshds ** (0.1D1 / 0.3D1)
      u9 = u8 ** 2
      u10 = 0.1D1 / u9
      contrf1DTD(i) =contrf1DTD(i) + (u6 * u7 * u10)*casp
      u12 = u3 * u5
      u22 = 0.1D1 / u9 / permshds
      ContrDOWN1dTD(j,i)= ContrDOWN1dTD(j,i) + (0.5D1 / 0.3D1 * u12 * (HauteurMan * dLmoyenDOWN + dHauteurManDOWN * Lmoyendiag) * u7 * u10 - 0.2D1 / 0.3D1 * u6 * dpermshdsDOWN * u7 * u22)*casp
      ContrUP1dTD(j,i) = ContrUP1dTD(j,i)+ (0.5D1 / 0.3D1 * u12 * (HauteurMan * dLmoyenUP + dHauteurManUP * Lmoyendiag) * u7 * u10 - 0.2D1 / 0.3D1 * u6 * dpermshdsUP * u7 * u22)*casp
!      if (temps.gt.63.75d0) then
!print*, "hauteurman", hauteurman
!print*, "nmanlocal", nmanlocal
!print*, "pentelocal", pentelocal
!print*, "dlp",dlp
!print*, "dpm", dpm
!print*, "lmoyendiag", lmoyendiag
!print*, "permshds", permshds
!print*,"dpermShdsdown",dpermShdsdown
!print*, "dpermshdsup", dpermshdsup
! print*, "dLmoyenUP", dLmoyenUP
! print*,  "dLmoyenDOWN", dLmoyenDOWN
!
!print*, "result",(u6 * u7 * u10)*casp
!pause
!endif
    
    exit 
    endif
enddo



end

!-----------------------------------------------------------------------------------------------------
subroutine comp_nb_v1d
use maillage_1d, only : nn1d,nb_v1d_m,nsources,nconnect
use riviere, only :  noeudsource, grandeurimposee
 implicit none
 integer i,j,cal_ne
cal_ne=0.d0

do i=1,nn1d
    do j=1,nsources
        if ((i.eq.noeudsource(j)).and.(grandeurimposee(j).eq."H")) then
            cal_ne=cal_ne+1
            goto 452
        endif
    enddo
            cal_ne=cal_ne+1+nconnect(i)
 452 continue   
enddo

nb_v1d_m=cal_ne

 end
 
  !-------------------------------------------------------------------------------------------------------------------     
        
               subroutine GestionDT
        use Observation1d
        use newtonp, only:   NewtonIT,lowlimIT,AugmentDT,highlimit,reducdt2,dtmin,dtmax
        use temporal, only : temps, dt
        implicit none
        integer i25
        
                        !Recalibration de dt en fonction de la rapidité de la convergence
                if (NewtonIT.le.lowlimIT) Then
                        dt=dt*Augmentdt
                elseif  (NewtonIT.ge.highlimIT) Then
                        dt=dt*reducDT2
                endif                    
                if (dt<dtmin) then
                        dt=dtmin
                elseif (dt.ge.dtmax) then
                        dt=dtmax
                endif
                
                     !Recalibration pour atteindre les points d'observation
             if (nextobs.le.nobst) then               
                                !Recalcul de l'indice du prochain point d'observation
                                nextobs=nobst+1
                                do i25=nobst,1,-1
                                    if ((temps)<obst(i25)) then
                                    nextobs=nextobs-1
                                    else
                                    exit
                                    endif
                                enddo                               
    
                                if (nextobs.le.nobst) then      
                                        if ((temps+dt).ge.obst(nextobs)) then
                                        dt=obst(nextobs)-temps
                                        elseif ((temps+2*dt).ge. obst(nextobs)) then
                                        dt=(obst(nextobs)-temps)/2.0d0
                                        endif
                                endif
            endif   
           print*,""   
                      
             end subroutine
             
!-------------------------------------------------------------------------------------------------------------------     
               subroutine GestionDTNONCONV
               !Modifie dt en cas de non convergence en mmax iteration. Puis, recalibraiton eventuelle du dt en fonction des temps d'observation
        use temporal, only:   temps,dt
        use Newtonp, only:   dtmin,reducdt1,dtmin,reducdt1
        use resolution1d, only : omega
        use observation1d
        implicit none     
        
                                                                    temps=temps-dt
                                                                    dt=dmax1(dtmin,dt*reducdt1)                                                                   
                                                                    if (nextobs.le.nobst) then
                                                                        if (obst(nextobs).gt.0.d0) then
                                                                            if ((temps+2*dt).ge. obst(nextobs)) then
                                                                                    dt=(obst(nextobs)-temps)/2.0d0
                                                                            endif
                                                                        endif
                                                                    endif
                                                                    if (dt.le. dtmin) then
                                                                    omega=1.0d0
                                                                    endif  
                                                                    temps=temps+dt                 

             end subroutine
!-------------------------------------------------------------------------------------------------------------------     
!-------------------------------------------------------------------------------------------------------------------     
              
             Subroutine NDMAIN (AX,AI,b,x,n,ne,iter,iter_tot)
             
!
   INTEGER NMAX, NEMAX, LVALUE, LINDEX,iter_tot
   PARAMETER  (NMAX=10000,NEMAX=20000,LVALUE=10000000,LINDEX=1000000)
   INTEGER KEEP_F (20), INDEXE_F (LINDEX), INFO (40),I, ICNTL (20), N, NE,AI (*)   !AI (2*NEMAX)
!
   DOUBLE PRECISION B (*), X (*), W (4*NMAX), VALUE (LVALUE), AX (*)  !AX (NEMAX)
   DOUBLE PRECISION CNTL (10), RINFO (20)
    
			!----------------que la résolution----------la matrice n'a pas changé---------------------!
  IF(iter>1)  then
     CALL UMD2SO (N, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEXE_F,KEEP_F, B, X, W, CNTL, ICNTL, INFO, RINFO)
     IF (INFO (1) .LT. 0) STOP
     return
   endif
                        !-----------------------------------------------------------------------------------------!
   DO I = 1, NE
     INDEXE_F(I) = AI (I)
     INDEXE_F(NE+I) = AI (NE+I)
     VALUE (I) = AX (I)
   ENDDO
			!---------------- la résolution----------la matrice a completement changé-----------------!
   IF(iter_tot==1) then
     CALL UMD21I (KEEP_F, CNTL, ICNTL)
     icntl(1)=61
     icntl(2)=61
     ICNTL (3) = 2
     icntl(6)=1
     CALL UMD2FA (N, NE, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEXE_F,KEEP_F, CNTL, ICNTL, INFO, RINFO)
  else
			!---------------- la résolution----------la matrice a changé uniquement les valeurs-------!
      CALL UMD2RF (N, NE, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEXE_F,KEEP_F, CNTL, ICNTL, INFO, RINFO)
  endif

   IF (INFO (1) .LT. 0) STOP
   ICNTL (3) = 2
     CALL UMD2SO (N, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEXE_F,KEEP_F, B, X, W, CNTL, ICNTL, INFO, RINFO)
   IF (INFO (1) .LT. 0) STOP

   return
   END
   
!-------------------------------------------------------------------------------------------------------------------     
!-------------------------------------------------------------------------------------------------------------------     
subroutine calc_flux1d
 use maillage_1D,only : nsources,nconnect,nconnected,mconnected,lengthM_1D,lengthM_1Dx,mn1D,nm1D,nn1d,cn
 use ecoulement1d, only : sommepluie1d,qpsrain1d,flux1d,widthmoyen_Maille,rh,dwidthprime,hauteur,amont, aval,kp,dkp,widthPRIME,widthmoyen,HplusZ,Number_inflow,Number_outflow,noeud_inflow,maille_inflow,noeud_outflow,maille_outflow,grad_HplusZ
 use riviere, only : condkp,hds1d, noeudsource,width,nman,pentel,penter,grandeurimposee,valeurimposee
 use temporal, only : temps,dt
 use observation1d, only : noeud_obs,nbr_noeuds_obs
 use bilanmasse1d, only : bilanq1d,bilanQpsRAIN1d, bilanqSourceNEUMANN
 use couplage, only : Coef_Asperites1d_1
 implicit none
 integer i,j,noeudamont,noeudaval,o79,o80,mailleconnec,i25
 real*8 lengthloc,penteldown,pentelup,penterup,penterdown,zup,zdown,ldown,lup,hauteurup,hauteurdown,lgthmail,nmanlocal,KPlo

     
do i=1,nn1d
      !H plus Z
     HplusZ(i)=cn(3,i)+hauteur(i) 
       !Lprime
     widthprime(i)=width(i)  +  (dmax1(hauteur(i)-hds1d,0.d0))* dwidthprime(i)
     !Lmoyen
     widthmoyen(i)= 0.5d0*(  width(i)+widthprime(i))
      !Rayon Hydraulique
     Rh(i)=((dmax1(hauteur(i)-hds1d,0.d0))*widthmoyen(i) )/  ( width(i) +( dmax1(hauteur(i)-hds1d,0.d0))*( 1.0d0/dSIN(pentel(i))+1.0d0/dSIN(penteR(i)))  )
enddo
  !Moyenne sur les deux noeuds d'une maille de Lmoyen 
 ! print*, nm1d
 do i=1,nm1d
 !print*, i
    widthmoyen_Maille(i)= 0.5d0*(  widthmoyen(mn1D(1,i))+widthmoyen(mn1D(2,i))  )
 enddo  
 !Détermination du noeud amont pour chaque maille
 do i=1,nm1d
     if (HplusZ(mn1D(1,i)).gt.HplusZ(mn1D(2,i))) then
         amont(i)=mn1D(1,i)
         aval(i)=mn1D(2,i)
     else
         amont(i)=mn1D(2,i) 
         aval(i)=mn1D(1,i)
     endif
 enddo
 
!Gradient de H+Z
 do i=1,nn1d
     do j=1,nconnect(i)
     mailleconnec=mconnected(j,i)
     lengthloc=lengthM_1D(mailleconnec)
        grad_HplusZ(j,i)=  -1.0d0*dabs((-HplusZ(i)+HplusZ(nconnected(j,i)))/lengthloc)
     enddo
 enddo

 
 !Calcul des flux1d et des flux1d cumulés pour tous les noeuds intérieurs, et aussi pour les sources a hauteur imposee
  do i=1,nn1d
     do j=1,nconnect(i)     
            noeudamont=amont(mconnected(j,i))
            noeudaval=aval(mconnected(j,i)) 
!            if (i.eq.1) then
!            print*, "am", noeudamont,"hauteur" ,hauteur(noeudamont)
!            print*, "av", noeudaval,"hauteur" ,hauteur(noeudaval)
!            endif
            Nmanlocal=(Nman(noeudamont)+Nman(noeudaval))*0.5d0  
            KPlo=(widthmoyen(noeudamont)*dmax1(hauteur(noeudamont)-hds1d,0.d0)*(rh(noeudamont))**(2.0d0/3.0d0) ) /(nmanlocal*dsqrt(dabs( grad_HplusZ(j,i)  )))  
!            if (i.eq.1) then
!            print*, "kp", kplo
!            print*, "casp", Coef_Asperites1d_1 (j,i)
!            endif          
            if ((KPlo.gt.condkp).or.(hauteur(noeudamont).le.hds1d)) then
                flux1d(j,i)=0.d0
            else
                kplo=kplo* Coef_Asperites1d_1 (j,i)
                flux1d(j,i)=  -1.0d0*kplo*grad_HplusZ(j,i)
            endif
            if (amont(mconnected(j,i)).eq. i)then
                 bilanq1d(i)=bilanq1d(i)-dt*flux1d(j,i)
            else if  (amont(mconnected(j,i)).eq. nconnected(j,i)) then
                 bilanq1d(i)=bilanq1d(i)+dt*flux1d(j,i)
            else
                print*, "error flux1d"
                pause
                stop
            endif
       enddo
       !print*, i, bilanq1d(i)
 enddo
 
  !bilan sur la source considérée comme une condition de Neumann
 if (nsources.gt.0) then
     do i=1,nsources
     if (grandeurimposee(i).eq."F") then
        bilanqSourceNEUMANN(noeudsource(i))=bilanqSourceNEUMANN(noeudsource(i))+valeurimposee(i)*dt
     endif
     enddo
 endif
 
 do i=1,nn1d    
     widthprime(i)=width(i)  +  ((hauteur(i))* dwidthprime(i) )
 enddo 
 !bilan sur qpsrain1d
 !Bilan des tes termes puits sources
 do i=1,nm1d
 noeudamont=amont(i)
 noeudaval=aval(i)
 bilanQpsRAIN1d(i)=bilanQpsRAIN1d(i)+lengthM_1D(i)*qpsrain1d(i)*dt*(0.5d0*widthprime(mn1d(1,i))+ 0.5d0*widthprime(mn1d(2,i))) 
 do i25=1,nsources
            if ((noeudamont.eq.noeudsource(i25).and.grandeurimposee(i25).eq."H").or.(noeudaval.eq.noeudsource(i25).and.grandeurimposee(i25).eq."H")) then
             bilanQpsRAIN1d(i)=bilanQpsRAIN1d(i)-0.5d0*lengthM_1D(i)*qpsrain1d(i)*dt*(0.5d0*widthprime(mn1d(1,i))+ 0.5d0*widthprime(mn1d(2,i))) 
             exit
            endif
enddo
enddo
 


 !ECRITURE DES flux1d
 !Ecriture de la première ligne au temps t=0
 if (temps.eq.dt) then
  write(1950,'(a10)',advance='no')"temps"
 do o79=1,nbr_noeuds_obs
  write(1950,'(a10)',advance='no') "noeudOBS"
  write(1950,'(a12)',advance='no') "NbrVoisins"
 do o80=1,nconnect(noeud_obs(o79))
 write(1950,'(a10)',advance='no') "in/out"
 write(1950,'(a10)',advance='no') "flux1d"
 enddo
 enddo
  write(1950,*) ""
 endif
 !Ecriture du flux1d à tous les pas de temps
 write(1950,'(f20.10)',advance='no')temps
 do o79=1,nbr_noeuds_obs
 if (o79.ne.nbr_noeuds_obs) then
     write(1950,'(i8,i8)',advance='no')noeud_obs(o79),nconnect(noeud_obs(o79))
     !print*, "noeudobs", noeud_obs(o79)
    ! print*, "noeuds voisins", nconnected(:,noeud_obs(o79))
    ! print*, "nconnect",nconnect(noeud_obs(o79))
     do o80=1,nconnect(noeud_obs(o79))
    ! print*, "o80",o80
    ! print*,mconnected(o80,noeud_obs(o79))
    ! print*, "amont",amont(mconnected(o80,noeud_obs(o79)))
     if (amont(mconnected(o80,noeud_obs(o79))).eq. noeud_obs(o79))then
     write(1950,'(a10,f20.10)',advance='no') "outflow",flux1d(o80,noeud_obs(o79))
     else if  (amont(mconnected(o80,noeud_obs(o79))).eq. nconnected(o80,noeud_obs(o79))) then
     write(1950,'(a10,f20.10)',advance='no') "inflow",flux1d(o80,noeud_obs(o79))
     else
     print*, "error flux1d"
     pause
     stop
     endif
     enddo
else
     write(1950,'(i8,i8)',advance='no')noeud_obs(o79),nconnect(noeud_obs(o79))
     !print*, "noeudobs", noeud_obs(o79)
    ! print*, "noeuds voisins", nconnected(:,noeud_obs(o79))
    ! print*, "nconnect",nconnect(noeud_obs(o79))
     do o80=1,nconnect(noeud_obs(o79))
     if (o80.ne.nconnect(noeud_obs(o79))) then
        ! print*, "o80",o80
        ! print*,mconnected(o80,noeud_obs(o79))
        ! print*, "amont",amont(mconnected(o80,noeud_obs(o79)))
         if (amont(mconnected(o80,noeud_obs(o79))).eq. noeud_obs(o79))then
         write(1950,'(a10,f20.10)',advance='no') "outflow",flux1d(o80,noeud_obs(o79))
         else if  (amont(mconnected(o80,noeud_obs(o79))).eq. nconnected(o80,noeud_obs(o79))) then
         write(1950,'(a10,f20.10)',advance='no') "inflow",flux1d(o80,noeud_obs(o79))
         else
         print*, "error flux1d"
         pause
         stop
         endif
    else
           ! print*, "o80",o80
        ! print*,mconnected(o80,noeud_obs(o79))
        ! print*, "amont",amont(mconnected(o80,noeud_obs(o79)))
         if (amont(mconnected(o80,noeud_obs(o79))).eq. noeud_obs(o79))then
         write(1950,'(a10,f20.10)') "outflow",flux1d(o80,noeud_obs(o79))
         else if  (amont(mconnected(o80,noeud_obs(o79))).eq. nconnected(o80,noeud_obs(o79))) then
         write(1950,'(a10,f20.10)') "inflow",flux1d(o80,noeud_obs(o79))
         else
         print*, "error flux1d"
         pause
         stop
         endif
        endif
     enddo
     endif
 enddo
  !pause
 end subroutine
 
!-------------------------------------------------------------------------------------------------------------------      
!-------------------------------------------------------------------------------------------------------------------       
 subroutine Plotting
 use maillage_1d, only : cn,nn1d,mn1d
 use ecoulement1d, only : hauteur
 use observation1d, only :nobst,obst,noeud_obs,nbr_noeuds_obs,nbr_mailles_cannal,numero_maille_cannal,ncannaux,canal,lengthcanalamont
 use temporal, only : temps
 implicit none
 
 
 integer i78,i52,i53,i54
 
 !ECRITURE DES HAUTEURS A CHAQUE PAS DE TEMPS pour les point considérés
 if (temps.eq.0.d0) then
 do i78=1,nbr_noeuds_obs
 if (i78.ne.nbr_noeuds_obs) then
 write(191,'(a20,i20)',advance='no') "temps",noeud_obs(i78)
 else
  write(191,'(a20,i20)') "temps",noeud_obs(i78)
  endif
 enddo
 endif 
 do i78=1,nbr_noeuds_obs
 if (i78.ne.nbr_noeuds_obs) then
    write(191,'(2f20.10)',advance='no')temps, hauteur(noeud_obs(i78))
    else
    write(191,'(2f20.10)')temps, hauteur(noeud_obs(i78))
    endif
 enddo
 
 ! ECRITURE du profil général de hauteurs, uniquement aux temps considérés
                do i78=1,nobst                
                        if ((temps.eq.obst(i78))) then
                                call hauteurs_finales
                        
                                write(192,'(a10,f20.10)')"temps",temps
                                write(192,*)"noeud altitude cannal lengthAMONT Hauteur"
                                do i53=1,ncannaux
                                    do i52=1,nn1d                                    
                                         if (canal(i52).eq.i53) then
                                             write(192,'(i8,f20.10,i8,2f20.10)')i52,cn(3,i52),canal(i52),lengthcanalamont(i52),hauteur(i52)
                                         endif
                                    enddo      
                                enddo     
                                
                                
                                                           
                                exit                        
                        endif
                enddo               
end subroutine
!-------------------------------------------------------------------------------------------------------------------      
!------------------------------------------------------------------------------------------------------------------- 
 subroutine calc_volume1d
use maillage_1D,only : lengthM_1D,mn1D,nm1D,nn1d,cn
 use ecoulement1d, only : widthmoyen_Maille,dwidthprime,hauteur,widthPRIME,widthmoyen
 use riviere, only : width,pentel,penter
 use temporal, only : temps
 use bilanmasse1d, only : volume1d_INI,volume1d
 implicit none
 integer i
volume1d=0.d0

do i=1,nn1d
       !Lprime
     widthprime(i)=width(i)  + hauteur(i)* dwidthprime(i)
     !Lmoyen
     widthmoyen(i)= 0.5d0*(  width(i)+widthprime(i))
enddo
  !Moyenne sur les deux noeuds d'une maille de Lmoyen 
 do i=1,nm1d
 !print*, i
    widthmoyen_Maille(i)= 0.5d0*(  widthmoyen(mn1D(1,i))+widthmoyen(mn1D(2,i))  )
 enddo  
 
 do i=1,nm1d
 volume1d=volume1d+ widthmoyen_Maille(i)*(0.5d0*hauteur(mn1d(1,i))+0.5d0*hauteur(mn1d(2,i)))*lengthM_1D(i)
 enddo

 if (temps.eq.0.d0) then
 volume1d_ini=volume1d
 endif
 
end subroutine
!-------------------------------------------------------------------------------------------------------------------      
!------------------------------------------------------------------------------------------------------------------- 
 subroutine Bilan_masse1d
 use maillage_1d, only : nm1d,nsources,nexits,nn1d
 use riviere,only : noeudsource,noeudsortie,grandeurimposee
 use temporal, only : temps
 use bilanmasse1d
 implicit none
 real*8 integ_hauteur,integ_q_entree,integ_Q_sortie,arraybilan(3)
 integer i24,i25
INTEG_hauteur=0.d0
INTEG_Q_entree=0.d0
INTEG_Q_sortie=0.d0

!Delta volume1d
INTEG_hauteur= volume1d-volume1d_ini
!Ajout de la pluie1d
integ_q_entree=integ_q_entree+ sum(bilanQpsRAIN1d(1:nm1d))
!AJout des sources
do i24=1,nn1d
    if (nsources.gt.0) then
        do i25=1,nsources
            if (i24.eq.noeudsource(i25).and.grandeurimposee(i25).eq."H") then
                integ_q_entree=integ_q_entree-bilanq1d(i24)
!                print*, "toto"
!                print*, "i", i24
!                print*, bilanq1d(i24)
                goto 854
            elseif (i24.eq.noeudsource(i25).and.grandeurimposee(i25).eq."F") then
              integ_q_entree=integ_q_entree+bilanqSourceNEUMANN(i24)
              goto 854            
            endif  
        enddo  
    endif  
    
!Calcul des flux1d sortants
    if (nexits.gt.0) then
        do i25=1,nexits
            if (i24.eq.noeudsortie(i25)) then
                integ_q_sortie=integ_q_sortie+bilanq1d(i24)
                goto 854
            endif  
        enddo
    endif
    854 continue
enddo

!ordonner les 3 termes dans l'ordre croissant
arraybilan(1)=integ_hauteur
arraybilan(2)= integ_q_entree
arraybilan(3)= integ_q_sortie
!print*, arraybilan(1),"h"
!print*, arraybilan(2),"e"
!print*, arraybilan(3),"s"
!pause
call sort_REAL8array(arraybilan,3)

!Computer le bilan de masse final
!print*,INTEG_hauteur,INTEG_Q_entree,INTEG_Q_sortie
if(dabs(arraybilan(1))+dabs(arraybilan(2)).gt.dabs(arraybilan(3))) then
write(194,'(a15,2f20.10,a2)'),"temps",temps, 100.0d0*(((dabs(arraybilan(1))+dabs(arraybilan(2)))/dabs(arraybilan(3)))-1.0d0),"%"
elseif(dabs(arraybilan(1))+dabs(arraybilan(2)).le.dabs(arraybilan(3))) then
write(194,'(a15,2f20.10,a2)'),"temps",temps, 100.0d0*((dabs(arraybilan(3))/(dabs(arraybilan(1))+dabs(arraybilan(2))))-1.0d0),"%"
endif
end subroutine 

!------------------------------------------------------------------------------------------------------------------- 
subroutine interpolSOURCE
use maillage_1d,only : nsources,interpolerlasource,ntimezonessources
use riviere, only : slopesource,timeofzonesource,valeurimposee,valimposezone,noeudsource, grandeurimposee
use temporal, only : temps
implicit none
integer i31,j31

 if (nsources>0) then 
        if (interpolerlasource.eq..TRUE.) then
         do i31=1,ntimezonesSOURCES     
             if (i31==ntimezonesSOURCES) then
                     do j31=1,nsources
                        valeurimposee(j31)=valimposezone(j31,i31)                        
                     enddo
                      exit
             elseif ((temps >= timeofzonesource(i31)) .and.(temps < timeofzonesource(i31+1)))  then
                 do j31=1,nsources
                     slopesource(j31)=(valIMPOSEZone(j31,i31+1)-valIMPOSEZone(j31,i31))/(timeofzoneSOURCE(i31+1)-timeofzoneSOURCE(i31))
                     valeurimposee(j31)=valIMPOSEZone(j31,i31)+slopesource(j31)*(temps-timeofzoneSOURCE(i31))                                     
                 enddo
                 exit
             endif
         enddo 
        else
          do i31=1,ntimezonessources     
             if (i31==ntimezonessources) then
                     do j31=1,nsources
                        valeurimposee(j31)=valIMPOSEZone(j31,i31)                        
                      enddo
                      exit
             elseif ((temps >= timeofzoneSOURCE(i31)) .and.(temps .le. timeofzoneSOURCE(i31+1)))  then
                     do j31=1,nsources
                        valeurimposee(j31)=valIMPOSEZone(j31,i31)                       
                      enddo
                      exit
             endif
         enddo
        endif
endif
!pause
!print*,"caca"
!do i31=1,nsources
!print*, "temps",temps,"noeud",noeudsource(i31),"v", valeurimposee(i31)
!enddo
!pause
return
end subroutine

!------------------------------------------------------------------------------------------------------------------- 
subroutine sort_REAL8array(xx,taille)
IMPLICIT  NONE
      real*8 :: xx(*)
      real*8                                :: temp,minimum
      INTEGER                               :: taille
      INTEGER                               :: i,j
      INTEGER                               :: Location

      DO i = 1, taille-1			   
                 Minimum  = xx(i)
                 location=i
                 DO j = i+1, taille		
                     IF (xx(j) < Minimum) THEN	
                        Minimum  = xx(j)		!     
                        Location = j                
                     END IF
                 END DO
              temp= xx(i)
              xx(i)=xx(location)
              xx(location)=temp             
      END DO     
      return
   END SUBROUTINE  sort_REAL8array
!------------------------------------------------------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------------------------------- 
subroutine hauteurs_finales
use Initialisation1d, only : fichier
use maillage_1d, only : nn1d,lengthM_1DX,lengthM_1D,mconnected,cn,cosinusdepente2,nconnect,mconnected
use ecoulement1d, only: hauteur
use temporal, only : temps
IMPLICIT  NONE
integer i,j,ind
ind=index(fichier,' ')-1
         open(1951,file='Output/'//fichier(1:ind)//'_hauteurFINALES',status='replace')
write(1951,'(3a14,a19,f15.5)') "noeud","hauteur","charge", "temps = ", temps
do i=1,nn1d
    do j=1,nconnect(i)
    cosinusdepente2(j,i)=lengthM_1DX(mconnected(j,i))/lengthM_1D(mconnected(j,i))
    enddo  
enddo
     
do i=1,nn1d
write(1951,'(i8,2f20.10)')i, hauteur(i),cn(3,i)+hauteur(i)/((sum(cosinusdepente2(1:nconnect(i),i)))/nconnect(i))
enddo  
close (1951)
      return
   END SUBROUTINE  hauteurs_finales
   
   !------------------------------------------------------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------------------------------- 
subroutine Temps_simu
use autres, only : nb_periodes_sec,nb_periodes_max,nb_periodes_initial,nb_periodes_final,nb_periodes , temps_elapsed
use temporal, only: temps,start,finish
implicit none
if (temps.eq.0.d0) then
    call cpu_time(start) 
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, COUNT_MAX=nb_periodes_max)
    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
else
call cpu_time(finish)
CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
nb_periodes = nb_periodes_final - nb_periodes_initial
IF (nb_periodes_final < nb_periodes_initial) nb_periodes = nb_periodes + nb_periodes_max
temps_elapsed   = REAL(nb_periodes) / nb_periodes_sec
write(67,'(a10,f15.4,a9)')  "temps Réel" ,temps_elapsed, "secondes"
write(67,'(a10,f15.4,a9)')  "temps CPU " ,finish-start, "secondes"
endif    
return
end subroutine