
      implicit double precision (a-h,o-z)
      parameter (gamal=7./5., gamar=7./5., nx=100,ny=10)
	
      double precision den(1:nx,1:ny),velx(1:nx,1:ny),vely(1:nx,1:ny)
	double precision pre(1:nx,1:ny),amass(1:nx,1:ny)
      double precision sound(1:nx,1:ny),tau(1:nx,1:ny)
      double precision enin(1:nx,1:ny),Ener(1:nx,1:ny)

	double precision vertx(0:nx,0:ny),verty(0:nx,0:ny)
      double precision Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	double precision Pstar(-1:nx+1,-1:ny+1,8)
      
	double precision enthropy0(1:nx),enthropy1(1:nx)
      write(*,*) 'Input T:'
      read(*,*) supt
      	
      call initial(den,velx,vely,pre,sound,tau,enin,Ener,
     &	       vertx,verty,amass,gamal,gamar,nx,ny) 
      j=5
	do i=1,nx
	   enthropy0(i)=pre(i,j)/den(i,j)**gamar
	enddo
 

      open(unit=10,file='Inital_Mesh.plt',status='unknown')
      write(10,*) 'title="contour"'
      write(10,*) 'variables="x","y","z"'
      write(10,*)  'zone  i=', nx+1, 'j=' ,ny+1, 'f=point'
      do j=0,ny
	  do i=0,nx
         write(10,102) vertx(i,j),verty(i,j),1.
        enddo 
      enddo
  102 format(3(F10.6, 1X))

      dt0=9999.
      t=0.	                          
	do while(t.lt.supt)

      Vstarx=99999.
	Vstary=99999.
	Pstar =99999.      
      call innodesolver(den,velx,vely,pre,sound,vertx,verty,nx,ny,
     &	              Vstarx,Vstary,Pstar)
             
      call getdt(sound,amass,tau,vertx,verty,Vstarx,Vstary,nx,ny,
     &	             dt0,dt)                  
	dt=0.0004
	dt0=dt  

      if(t+dt.gt.supt) dt=supt-t
        t=t+dt
        write(*,*) "T=",t,"dt=",dt
      

	call update(den,velx,vely,pre,sound,tau,enin,Ener,
     & vertx,verty,Vstarx,Vstary,Pstar,amass,nx,ny,gamal,gamar,dt)
	
	enddo       
       
      write(*,*) 'T=',t

      open(unit=100,file='Result.plt',status='unknown')
      write(100,*) 'title="contour"'
      write(100,*) 'variables="x","den","velx","vely","pre" "enin"'
	j=5
      do i=1,nx
        write(100,101)  (vertx(i-1,j)+vertx(i,j))/2.,
     &        den(i,j),velx(i,j),vely(i,j),pre(i,j),enin(i,j)
      enddo 
  101 format(10(F10.6, 1X))


      open(unit=10,file='Result_Mesh.plt',status='unknown')
      write(10,*) 'title="contour"'
      write(10,*) 'variables="x","y","z"'
      write(10,*)  'zone  i=', nx+1, 'j=' ,ny+1, 'f=point'
      do j=0,ny
	  do i=0,nx
         write(10,102) vertx(i,j),verty(i,j),1.
	
        enddo 
      enddo

      j=5
	do i=1,nx
	   enthropy1(i)=pre(i,j)/den(i,j)**gamar
	   print*, enthropy1(i)/ enthropy0(i)
	enddo

      
     
      end
c------------------------------------------------------------------------------------
c               initial
c------------------------------------------------------------------------------------
      subroutine initial(den,velx,vely,pre,sound,tau,enin,Ener,
     &	       vertx,verty,amass,gamal,gamar,nx,ny)

      implicit double precision (a-h,o-z) 

      double precision den(1:nx,1:ny),velx(1:nx,1:ny),vely(1:nx,1:ny)
	double precision pre(1:nx,1:ny),amass(1:nx,1:ny)
      double precision sound(1:nx,1:ny),tau(1:nx,1:ny)
      double precision enin(1:nx,1:ny),Ener(1:nx,1:ny)
      double precision area(1:nx,1:ny)

	double precision vertx(0:nx,0:ny),verty(0:nx,0:ny)


      xl=0.
      xr=1.
      yd=0.
	yu=0.1
 
      dx=(xr-xl)/nx      !  空间步长
	dy=(yu-yd)/ny
    

      do  i=1,nx                ! 单元内信息的定义
	  do  j=1,ny	    
           if(i<=nx/2) then
              den(i,j)=1.
              velx(i,j)=0.
	        vely(i,j)=0.
              pre(i,j)=1.

           else
              den(i,j)=0.125
              velx(i,j)=0.
	        vely(i,j)=0.
              pre(i,j)=0.1
		 endif
      enddo
	enddo 
                                ! 节点坐标
      do i=0,nx
	   do j=0,ny
	      vertx(i,j)=i*dx
	      verty(i,j)=j*dy
	   enddo
	enddo

                               ! 单元的面积
	do i=1,nx
	   do j=1,ny
            area(i,j)=0.0

			    xr1=vertx(i-1,j-1)           
	            yr1=verty(i-1,j-1)
	            xr2=vertx(i,j-1)           
	            yr2=verty(i,j-1)	            
	         area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

			    xr1=vertx(i,j-1)           
	            yr1=verty(i,j-1)
	            xr2=vertx(i,j)           
	            yr2=verty(i,j)	            
	         area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

			    xr1=vertx(i,j)           
	            yr1=verty(i,j)
	            xr2=vertx(i-1,j)           
	            yr2=verty(i-1,j)	            
	         area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)

			    xr1=vertx(i-1,j)           
	            yr1=verty(i-1,j)
	            xr2=vertx(i-1,j-1)           
	            yr2=verty(i-1,j-1)	            
	         area(i,j)=area(i,j)+0.5*(xr1*yr2-xr2*yr1)
		     
	         amass(i,j)=den(i,j)*area(i,j)       !每个单元的质量
	   enddo
	enddo

      do  i=1,nx                ! 单元内信息的定义
	  do  j=1,ny	    
           if(i<=nx/2) then

              velx(i,j)=den(i,j)*area(i,j)*velx(i,j)/amass(i,j)
	        vely(i,j)=den(i,j)*area(i,j)*vely(i,j)/amass(i,j)

	        sound(i,j)=sqrt( gamal*pre(i,j)/den(i,j) )
	        enin(i,j)=pre(i,j)/(gamal-1.)/den(i,j)
	        tau(i,j)=1./den(i,j)
              Ener(i,j)=den(i,j)*area(i,j)/amass(i,j)*
     &			(enin(i,j)+0.5*((velx(i,j))**2+(vely(i,j))**2))
           else
              velx(i,j)=den(i,j)*area(i,j)*velx(i,j)/amass(i,j)
	        vely(i,j)=den(i,j)*area(i,j)*vely(i,j)/amass(i,j)

	        sound(i,j)=sqrt( gamar*pre(i,j)/den(i,j) )
	        enin(i,j)=pre(i,j)/(gamar-1.)/den(i,j)
	        tau(i,j)=1./den(i,j)
              Ener(i,j)=den(i,j)*area(i,j)/amass(i,j)*
     &			(enin(i,j)+0.5*((velx(i,j))**2+(vely(i,j))**2))
		 endif
        enddo
	enddo 

      end


c------------------------------------------------------------------------------------
c                internodesolver
c------------------------------------------------------------------------------------
	subroutine innodesolver(ghostden,ghostvelx,ghostvely,ghostpre,
     & ghostsound,ghostvertx,ghostverty,nx,ny,Vstarx,Vstary,Pstar)
	implicit double precision (a-h,o-z)

      double precision ghostden(1:nx,1:ny),ghostvelx(1:nx,1:ny)
	double precision ghostvely(1:nx,1:ny)
	double precision ghostpre(1:nx,1:ny),ghostsound(1:nx,1:ny)
	double precision ghostvertx(0:nx,0:ny),ghostverty(0:nx,0:ny)

c      以下才是虚拟的网格名称
      double precision den(0:nx+1,0:ny+1),velx(0:nx+1,0:ny+1)
	double precision vely(0:nx+1,0:ny+1)
	double precision pre(0:nx+1,0:ny+1),sound(0:nx+1,0:ny+1)
	double precision vertx(-1:nx+1,-1:ny+1),verty(-1:nx+1,-1:ny+1)

      double precision Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	double precision Pstar(-1:nx+1,-1:ny+1,8)   ! 直接定义在单元内

!可调节参数
      z=1.
      
c     虚拟网格拓扑
	do i=0,nx
	   do j=0,ny
	      vertx(i,j)=ghostvertx(i,j)
	      verty(i,j)=ghostverty(i,j)
	   enddo
	enddo
      
	j=-1  !下边界
      do i=0,nx
	      vertx(i,j)=ghostvertx(i,j+2)
	      verty(i,j)=2*ghostverty(i,j+1)-ghostverty(i,j+2)
      enddo     

	j=ny+1
      do i=0,nx
	      vertx(i,j)=ghostvertx(i,j-2)
	      verty(i,j)=2*ghostverty(i,j-1)-ghostverty(i,j-2)
      enddo     

	i=-1
      do j=0,ny
	      vertx(i,j)=2*ghostvertx(i+1,j)-ghostvertx(i+2,j)
	      verty(i,j)=ghostverty(i+2,j)
      enddo     

	i=nx+1
      do j=0,ny
	      vertx(i,j)=2*ghostvertx(i-1,j)-ghostvertx(i-2,j)
	      verty(i,j)=ghostverty(i-2,j)
      enddo     
    
	i=-1
	j=-1
      vertx(i,j)=2*ghostvertx(i+1,j+1)-ghostvertx(i+2,j+2)
	verty(i,j)=2*ghostverty(i+1,j+1)-ghostverty(i+2,j+2)

	i=-1
	j=ny+1
      vertx(i,j)=2*ghostvertx(i+1,j-1)-ghostvertx(i+2,j-2)
	verty(i,j)=2*ghostverty(i+1,j-1)-ghostverty(i+2,j-2)

	i=nx+1
	j=-1
      vertx(i,j)=2*ghostvertx(i-1,j+1)-ghostvertx(i-2,j+2)
	verty(i,j)=2*ghostverty(i-1,j+1)-ghostverty(i-2,j+2)

	i=nx+1
	j=ny+1
      vertx(i,j)=2*ghostvertx(i-1,j-1)-ghostvertx(i-2,j-2)
	verty(i,j)=2*ghostverty(i-1,j-1)-ghostverty(i-2,j-2)

      
 ! 延拓虚拟网格中的量的拓扑
	do i=1,nx               
	   do j=1,ny
	      pre(i,j)=ghostpre(i,j)
	      den(i,j)=ghostden(i,j)
	      velx(i,j)=ghostvelx(i,j)
	      vely(i,j)=ghostvely(i,j)
	      sound(i,j)=ghostsound(i,j)
	   enddo
	enddo

      i=0                    ! 左边界   处理是否有问题
	do j=1,ny
	   pre(i,j)=ghostpre(i+1,j)
	   den(i,j)=ghostden(i+1,j)
	   velx(i,j)=-ghostvelx(i+1,j)
	   vely(i,j)=ghostvely(i+1,j)
	   sound(i,j)=ghostsound(i+1,j)
	enddo

      i=nx+1                ! 右边界
	do j=1,ny
	   pre(i,j)=ghostpre(i-1,j)
	   den(i,j)=ghostden(i-1,j)
	   velx(i,j)=-ghostvelx(i-1,j)
	   vely(i,j)=ghostvely(i-1,j)
	   sound(i,j)=ghostsound(i-1,j)
	enddo

      j=0                   ! 下边界
	do i=1,nx
	    pre(i,j)=ghostpre(i,j+1)
	    den(i,j)=ghostden(i,j+1)
	    velx(i,j)=ghostvelx(i,j+1)        
	    vely(i,j)=-ghostvely(i,j+1)
	    sound(i,j)=ghostsound(i,j+1)
	enddo

      j=ny+1                ! 上边界
	do i=1,nx
	    pre(i,j)=ghostpre(i,j-1)
	    den(i,j)=ghostden(i,j-1)
	    velx(i,j)=ghostvelx(i,j-1)
	    vely(i,j)=-ghostvely(i,j-1)
	    sound(i,j)=ghostsound(i,j-1)
	enddo
      
	i=0                   ! 左下顶点
	j=0
      pre(i,j)=ghostpre(i+1,j+1)
	den(i,j)=ghostden(i+1,j+1)
	velx(i,j)=-ghostvelx(i+1,j+1)         
	vely(i,j)=-ghostvely(i+1,j+1)
      sound(i,j)=ghostsound(i+1,j+1)
      
	
	i=nx+1                !右下顶点
	j=0
	    pre(i,j)=ghostpre(i-1,j+1)
	    den(i,j)=ghostden(i-1,j+1)
	    velx(i,j)=-ghostvelx(i-1,j+1)
	    vely(i,j)=-ghostvely(i-1,j+1)
	    sound(i,j)=ghostsound(i-1,j+1)
		      
	i=0                  ! 左上顶点
	j=ny+1
	    pre(i,j)=ghostpre(i+1,j-1)
	    den(i,j)=ghostden(i+1,j-1)
	    velx(i,j)=-ghostvelx(i+1,j-1)
	    vely(i,j)=-ghostvely(i+1,j-1)
	    sound(i,j)=ghostsound(i+1,j-1)
		     
	i=nx+1              ! 右上顶点
	j=ny+1
	    pre(i,j)=ghostpre(i-1,j-1)
	    den(i,j)=ghostden(i-1,j-1)
	    velx(i,j)=-ghostvelx(i-1,j-1)
	    vely(i,j)=-ghostvely(i-1,j-1)
	    sound(i,j)=ghostsound(i-1,j-1)
  

      do j=0,ny               ! 顶点扫描
	  do i=0,nx	   

		  	 temptA=0.0
	         temptB=0.0
	         temptC=0.0
	         SMx=0.0
	         SMy=0.0   
	    	         	         
c1               
	         x1=vertx(i,j)        ! 节点坐标
	         y1=verty(i,j)
	         x2=vertx(i+1,j)      ! 连接点坐标
	         y2=verty(i+1,j)
       
			 call Length(x1,y1,x2,y2,alength)            ! 边的长度               
			 call Normal(x1,y1,x2,y2,anormalx,anormaly)  ! 边的法向
	         
			 anormalxnn=-anormalx              ! 1796 第一行
	         anormalynn=-anormaly
               
               knus1i=i+1
	         knus1j=j
	         ki=i+1
	         kj=j+1	
	 		 
	 addalpha=den(knus1i,knus1j)*sound(knus1i,knus1j)*z                
     &                  +den(ki,kj)*sound(ki,kj)*z   

	         temptA=temptA+alength*addalpha*anormalxnn**2
	         temptB=temptB+alength*addalpha*anormalynn**2
	         temptC=temptC+alength*addalpha*anormalxnn*anormalynn

	         starvk=(  pre(knus1i,knus1j)-pre(ki,kj)
     &		+den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &      (velx(knus1i,knus1j)*anormalx+
     &               vely(knus1i,knus1j)*anormaly)
     &        +den(ki,kj)*sound(ki,kj)*z*
     &         (velx(ki,kj)*anormalx
     &             +vely(ki,kj)*anormaly) )/addalpha 
    
	         SMx=SMx+alength*addalpha*starvk*anormalx
               SMy=SMy+alength*addalpha*starvk*anormaly
c2	      
	         x1=vertx(i,j)        ! 节点坐标
	         y1=verty(i,j)
	         x2=vertx(i,j+1)      ! 连接点坐标
	         y2=verty(i,j+1)
       
			 call Length(x1,y1,x2,y2,alength)            ! 边的长度               
			 call Normal(x1,y1,x2,y2,anormalx,anormaly)  ! 边的法向
	         
			 anormalxnn=-anormalx              ! 1796 第一行
	         anormalynn=-anormaly
               
               knus1i=i+1
	         knus1j=j+1
	         ki=i
	         kj=j+1
			  
	 addalpha=den(knus1i,knus1j)*sound(knus1i,knus1j)*z                
     &                  +den(ki,kj)*sound(ki,kj)*z   

	         temptA=temptA+alength*addalpha*anormalxnn**2
	         temptB=temptB+alength*addalpha*anormalynn**2
	         temptC=temptC+alength*addalpha*anormalxnn*anormalynn

	         starvk=( pre(knus1i,knus1j)-pre(ki,kj)
     &		+den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &      (velx(knus1i,knus1j)*anormalx
     &                   +vely(knus1i,knus1j)*anormaly)
     &        +den(ki,kj)*sound(ki,kj)*z*
     &         (velx(ki,kj)*anormalx+
     &              vely(ki,kj)*anormaly) )/addalpha
    
	         SMx=SMx+alength*addalpha*starvk*anormalx
               SMy=SMy+alength*addalpha*starvk*anormaly
c3	      
	         x1=vertx(i,j)        ! 节点坐标
	         y1=verty(i,j)
	         x2=vertx(i-1,j)      ! 连接点坐标
	         y2=verty(i-1,j)
       
			 call Length(x1,y1,x2,y2,alength)            ! 边的长度               
			 call Normal(x1,y1,x2,y2,anormalx,anormaly)  ! 边的法向
	         
			 anormalxnn=-anormalx              ! 1796 第一行
	         anormalynn=-anormaly
               
               knus1i=i
	         knus1j=j+1
	         ki=i
	         kj=j
			  
	      addalpha=den(knus1i,knus1j)*sound(knus1i,knus1j)*z                
     &                  +den(ki,kj)*sound(ki,kj)*z   

	         temptA=temptA+alength*addalpha*anormalxnn**2
	         temptB=temptB+alength*addalpha*anormalynn**2
	         temptC=temptC+alength*addalpha*anormalxnn*anormalynn

	         starvk=( pre(knus1i,knus1j)-pre(ki,kj)
     &		+den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &      (velx(knus1i,knus1j)*anormalx+
     &           vely(knus1i,knus1j)*anormaly)
     &        +den(ki,kj)*sound(ki,kj)*z*
     &         (velx(ki,kj)*anormalx
     &           +vely(ki,kj)*anormaly) )/addalpha
    
	         SMx=SMx+alength*addalpha*starvk*anormalx
               SMy=SMy+alength*addalpha*starvk*anormaly
	      
c4
	         x1=vertx(i,j)        ! 节点坐标
	         y1=verty(i,j)
	         x2=vertx(i,j-1)      ! 连接点坐标
	         y2=verty(i,j-1)
       
			 call Length(x1,y1,x2,y2,alength)            ! 边的长度               
			 call Normal(x1,y1,x2,y2,anormalx,anormaly)  ! 边的法向
	         
			 anormalxnn=-anormalx              ! 1796 第一行
	         anormalynn=-anormaly
               
               knus1i=i
	         knus1j=j
	         ki=i+1
	         kj=j
			  
	      addalpha=den(knus1i,knus1j)*sound(knus1i,knus1j)*z                
     &                  +den(ki,kj)*sound(ki,kj)*z   

	         temptA=temptA+alength*addalpha*anormalxnn**2
	         temptB=temptB+alength*addalpha*anormalynn**2
	         temptC=temptC+alength*addalpha*anormalxnn*anormalynn

	         starvk=(pre(knus1i,knus1j)-pre(ki,kj)
     &		+den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &      (velx(knus1i,knus1j)*anormalx+
     &              vely(knus1i,knus1j)*anormaly)
     &        +den(ki,kj)*sound(ki,kj)*z*
     &         (velx(ki,kj)*anormalx
     &               +vely(ki,kj)*anormaly) )/addalpha
    
	         SMx=SMx+alength*addalpha*starvk*anormalx
               SMy=SMy+alength*addalpha*starvk*anormaly

	      Vstarx(i,j)=(SMx*temptB-SMy*temptC)
     &		          /(temptA*temptB-temptC**2)
         	  Vstary(i,j)=(SMy*temptA-SMx*temptC)
     &		          /(temptA*temptB-temptC**2)  
c1	      
               x1=vertx(i,j)
	         y1=verty(i,j)
               x2=vertx(i+1,j)
	         y2=verty(i+1,j)
               knus1i=i+1        ! 公共边的上单元
	         knus1j=j
	         ki=i+1            ! 公共边的上单元
	         kj=j+1
	         
	         call Normal(x1,y1,x2,y2,anormalx,anormaly)

                Pstar(ki,kj,1)=pre(ki,kj)
     &			  +den(ki,kj)*sound(ki,kj)*z* 
     &         ((Vstarx(i,j)-velx(ki,kj))*anormalx +
     &          (Vstary(i,j)-vely(ki,kj))*anormaly )

	         Pstar(knus1i,knus1j,6)=pre(knus1i,knus1j)-      
     &	      den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &         ((Vstarx(i,j)-velx(knus1i,knus1j))*anormalx +
     &         (Vstary(i,j)-vely(knus1i,knus1j))*anormaly)
               p1=Pstar(ki,kj,1)
			 p6=Pstar(knus1i,knus1j,6)

               x1=vertx(i,j)
	         y1=verty(i,j)
               x2=vertx(i,j+1)
	         y2=verty(i,j+1)
               knus1i=i+1        ! 公共边的下单元
	         knus1j=j+1
	         ki=i            ! 公共边的上单元
	         kj=j+1
	         
	         call Normal(x1,y1,x2,y2,anormalx,anormaly)

                Pstar(ki,kj,3)=pre(ki,kj)
     &			  +den(ki,kj)*sound(ki,kj)*z* 
     &         ((Vstarx(i,j)-velx(ki,kj))*anormalx +
     &          (Vstary(i,j)-vely(ki,kj))*anormaly )

	         Pstar(knus1i,knus1j,8)=pre(knus1i,knus1j)-      
     &	       den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &         ((Vstarx(i,j)-velx(knus1i,knus1j))*anormalx +
     &         (Vstary(i,j)-vely(knus1i,knus1j))*anormaly)
               p3=Pstar(ki,kj,3)
			 p8=Pstar(knus1i,knus1j,8)


               x1=vertx(i,j)
	         y1=verty(i,j)
               x2=vertx(i-1,j)
	         y2=verty(i-1,j)
               knus1i=i        ! 公共边的下单元
	         knus1j=j+1
	         ki=i            ! 公共边的上单元
	         kj=j
	         
	         call Normal(x1,y1,x2,y2,anormalx,anormaly)

                Pstar(ki,kj,5)=pre(ki,kj)
     &			  +den(ki,kj)*sound(ki,kj)*z* 
     &         ((Vstarx(i,j)-velx(ki,kj))*anormalx +
     &          (Vstary(i,j)-vely(ki,kj))*anormaly )

	         Pstar(knus1i,knus1j,2)=pre(knus1i,knus1j)-       
     &	      den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &         ((Vstarx(i,j)-velx(knus1i,knus1j))*anormalx +
     &         (Vstary(i,j)-vely(knus1i,knus1j))*anormaly)
               p5=Pstar(ki,kj,5)
			 p2=Pstar(knus1i,knus1j,2)
      
	
	         x1=vertx(i,j)
	         y1=verty(i,j)
               x2=vertx(i,j-1)
	         y2=verty(i,j-1)
               knus1i=i        ! 公共边的下单元
	         knus1j=j
	         ki=i+1            ! 公共边的上单元
	         kj=j
	         
	         call Normal(x1,y1,x2,y2,anormalx,anormaly)

                Pstar(ki,kj,7)=pre(ki,kj)
     &			  +den(ki,kj)*sound(ki,kj)*z* 
     &         ((Vstarx(i,j)-velx(ki,kj))*anormalx +
     &          (Vstary(i,j)-vely(ki,kj))*anormaly )

	         Pstar(knus1i,knus1j,4)=pre(knus1i,knus1j)-       
     &	       den(knus1i,knus1j)*sound(knus1i,knus1j)*z*
     &         ((Vstarx(i,j)-velx(knus1i,knus1j))*anormalx +
     &         (Vstary(i,j)-vely(knus1i,knus1j))*anormaly)
               p7=Pstar(ki,kj,7)
			 p4=Pstar(knus1i,knus1j,4)

c      print*,p8,p3
c      print*,p2,p5
c      print*,p4,p7
c      print*,p6,p1
c      pause 615       
	  enddo	
	enddo 

!      do i=45,55
!	  do j=5,5
!      print*,Pstar(i,j,1),Pstar(i,j,2)       
!      print*,Pstar(i,j-1,6),Pstar(i,j-1,5)       
!      print*,Pstar(i,j,3),Pstar(i,j,4)       
!      print*,Pstar(i+1,j,8),Pstar(i+1,j,7)       
!      print*,Pstar(i,j,5),Pstar(i,j,6)       
!      print*,Pstar(i,j+1,2),Pstar(i,j+1,1)       
!      print*,Pstar(i,j,7),Pstar(i,j,8)       
!      print*,Pstar(i-1,j,4),Pstar(i-1,j,3) 
   
!	enddo
!	enddo
!	pause 627 

	end



c------------------------------------------------------------------------------------
c                update
c------------------------------------------------------------------------------------
      subroutine update(den,velx,vely,pre,sound,tau,enin,Ener,
     & vertx,verty,Vstarx,Vstary,Pstar,amass,nx,ny,gammal,gammar,dt)
      implicit double precision (a-h,o-z)
	
      double precision den(1:nx,1:ny),velx(1:nx,1:ny),vely(1:nx,1:ny)
	double precision pre(1:nx,1:ny),amass(1:nx,1:ny)
      double precision sound(1:nx,1:ny),tau(1:nx,1:ny)
      double precision enin(1:nx,1:ny),Ener(1:nx,1:ny)

	double precision vertx(0:nx,0:ny),verty(0:nx,0:ny)
      double precision Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
	double precision Pstar(-1:nx+1,-1:ny+1,8)

      double precision tempttau(1:nx,1:ny),temptvelx(1:nx,1:ny)
      double precision temptEner(1:nx,1:ny),temptvely(1:nx,1:ny)

	do i=1,nx
	   do j=1,ny

            Sum1=0.
            Sum2x=0.
	      Sum2y=0.
	      Sum3=0.
c1
			    ir=i-1
	            jr=j-1
	            ir1=i         ! 第一个点
	            jr1=j-1
      
	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            xrn1=xr+dt*Vstarx(ir,jr)
	            yrn1=yr+dt*Vstary(ir,jr)
	            xr1n1=xr1+dt*Vstarx(ir1,jr1)
	            yr1n1=yr1+dt*Vstary(ir1,jr1)

	            call Length(xr,yr,xr1,yr1,alengthn)
	            call Length(xrn1,yrn1,xr1n1,yr1n1,alengthn1)

                  call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
	       call Normal(xrn1,yrn1,xr1n1,yr1n1,anormalxn1,anormalyn1)
	            
	anormalxn=-anormalxn          ! 外法向
      anormalyn=-anormalyn
	anormalxn1=-anormalxn1
	anormalyn1=-anormalyn1
                  
                  Sum1=Sum1+( alengthn*anormalxn+alengthn1*anormalxn1 )*   
     &                      ( Vstarx(ir,jr)+Vstarx(ir1,jr1) ) +	
     &                      (alengthn*anormalyn+alengthn1*anormalyn1)*   
     &                      ( Vstary(ir,jr)+Vstary(ir1,jr1) ) 
      
	            Sum2x=Sum2x+alengthn*anormalxn*
     &				(Pstar(i,j,1)+Pstar(i,j,2))
				Sum2y=Sum2y+alengthn*anormalyn*
     &				(Pstar(i,j,1)+Pstar(i,j,2))
				  
				Sum3=Sum3+alengthn*(  Pstar(i,j,1)*
     & (Vstarx(ir,jr)*anormalxn+Vstary(ir,jr)*anormalyn)   + 			          
     &  Pstar(i,j,2)*
     & (Vstarx(ir1,jr1)*anormalxn+Vstary(ir1,jr1)*anormalyn)  )

c2
                  ir=i
	            jr=j-1
	            ir1=i         ! 第一个点
	            jr1=j
      
	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            xrn1=xr+dt*Vstarx(ir,jr)
	            yrn1=yr+dt*Vstary(ir,jr)
	            xr1n1=xr1+dt*Vstarx(ir1,jr1)
	            yr1n1=yr1+dt*Vstary(ir1,jr1)

	            call Length(xr,yr,xr1,yr1,alengthn)
	            call Length(xrn1,yrn1,xr1n1,yr1n1,alengthn1)

                  call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
	       call Normal(xrn1,yrn1,xr1n1,yr1n1,anormalxn1,anormalyn1)
	            
	anormalxn=-anormalxn          ! 外法向
      anormalyn=-anormalyn
	anormalxn1=-anormalxn1
	anormalyn1=-anormalyn1
                  
                  Sum1=Sum1+( alengthn*anormalxn+alengthn1*anormalxn1 )*   
     &                      ( Vstarx(ir,jr)+Vstarx(ir1,jr1) ) +	
     &                      (alengthn*anormalyn+alengthn1*anormalyn1)*   
     &                      ( Vstary(ir,jr)+Vstary(ir1,jr1) ) 
      
	            Sum2x=Sum2x+alengthn*anormalxn*
     &				(Pstar(i,j,3)+Pstar(i,j,4))
				Sum2y=Sum2y+alengthn*anormalyn*
     &				(Pstar(i,j,3)+Pstar(i,j,4))
				  
				Sum3=Sum3+alengthn*(  Pstar(i,j,3)*
     & (Vstarx(ir,jr)*anormalxn+Vstary(ir,jr)*anormalyn)   + 			          
     &  Pstar(i,j,4)*
     & (Vstarx(ir1,jr1)*anormalxn+Vstary(ir1,jr1)*anormalyn)  )

c3
	            ir=i
	            jr=j
	            ir1=i-1         ! 第一个点
	            jr1=j
      
	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            xrn1=xr+dt*Vstarx(ir,jr)
	            yrn1=yr+dt*Vstary(ir,jr)
	            xr1n1=xr1+dt*Vstarx(ir1,jr1)
	            yr1n1=yr1+dt*Vstary(ir1,jr1)

	            call Length(xr,yr,xr1,yr1,alengthn)
	            call Length(xrn1,yrn1,xr1n1,yr1n1,alengthn1)

                  call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
	       call Normal(xrn1,yrn1,xr1n1,yr1n1,anormalxn1,anormalyn1)
	            
	anormalxn=-anormalxn          ! 外法向
      anormalyn=-anormalyn
	anormalxn1=-anormalxn1
	anormalyn1=-anormalyn1
                  
                  Sum1=Sum1+( alengthn*anormalxn+alengthn1*anormalxn1 )*   
     &                      ( Vstarx(ir,jr)+Vstarx(ir1,jr1) ) +	
     &                      (alengthn*anormalyn+alengthn1*anormalyn1)*   
     &                      ( Vstary(ir,jr)+Vstary(ir1,jr1) ) 
      
	            Sum2x=Sum2x+alengthn*anormalxn*
     &				(Pstar(i,j,5)+Pstar(i,j,6))
				Sum2y=Sum2y+alengthn*anormalyn*
     &				(Pstar(i,j,5)+Pstar(i,j,6))
				  
				Sum3=Sum3+alengthn*(  Pstar(i,j,5)*
     & (Vstarx(ir,jr)*anormalxn+Vstary(ir,jr)*anormalyn)   + 			          
     &  Pstar(i,j,6)*
     & (Vstarx(ir1,jr1)*anormalxn+Vstary(ir1,jr1)*anormalyn)  )

c4
	            ir=i-1
	            jr=j
	            ir1=i-1         
	            jr1=j-1
      
	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            xrn1=xr+dt*Vstarx(ir,jr)
	            yrn1=yr+dt*Vstary(ir,jr)
	            xr1n1=xr1+dt*Vstarx(ir1,jr1)
	            yr1n1=yr1+dt*Vstary(ir1,jr1)

	            call Length(xr,yr,xr1,yr1,alengthn)
	            call Length(xrn1,yrn1,xr1n1,yr1n1,alengthn1)

                  call Normal(xr,yr,xr1,yr1,anormalxn,anormalyn)
	       call Normal(xrn1,yrn1,xr1n1,yr1n1,anormalxn1,anormalyn1)
	            
	anormalxn=-anormalxn          ! 外法向
      anormalyn=-anormalyn
	anormalxn1=-anormalxn1
	anormalyn1=-anormalyn1
                  
                  Sum1=Sum1+( alengthn*anormalxn+alengthn1*anormalxn1 )*   
     &                      ( Vstarx(ir,jr)+Vstarx(ir1,jr1) ) +	
     &                      (alengthn*anormalyn+alengthn1*anormalyn1)*   
     &                      ( Vstary(ir,jr)+Vstary(ir1,jr1) ) 
      
	            Sum2x=Sum2x+alengthn*anormalxn*
     &				(Pstar(i,j,7)+Pstar(i,j,8))
				Sum2y=Sum2y+alengthn*anormalyn*
     &				(Pstar(i,j,7)+Pstar(i,j,8))
				  
				Sum3=Sum3+alengthn*(  Pstar(i,j,7)*
     & (Vstarx(ir,jr)*anormalxn+Vstary(ir,jr)*anormalyn)   + 			          
     &  Pstar(i,j,8)*
     & (Vstarx(ir1,jr1)*anormalxn+Vstary(ir1,jr1)*anormalyn)  )
     	         
		    ! 求和运算完毕
                 tempttau(i,j)=tau(i,j)+Sum1*dt/4./amass(i,j)   ! 5.5
                 temptvelx(i,j)=velx(i,j)-Sum2x*dt/2./amass(i,j)  ! 5.6
	           temptvely(i,j)=vely(i,j)-Sum2y*dt/2./amass(i,j)
	           temptEner(i,j)=Ener(i,j)-Sum3*dt/2./amass(i,j)     ! 5.7
	   enddo
	enddo

      do i=0,nx
	   do j=0,ny
            vertx(i,j)=vertx(i,j)+dt*Vstarx(i,j)
            verty(i,j)=verty(i,j)+dt*Vstary(i,j)
	   enddo
	enddo
     
      do j=1,ny
      do i=1,nx
	   
	      tau(i,j)=tempttau(i,j)
	      velx(i,j)=temptvelx(i,j)
		  vely(i,j)=temptvely(i,j)
		  Ener(i,j)=temptEner(i,j) 
		  den(i,j)=1./tau(i,j)
            enin(i,j)=Ener(i,j)-0.5*(velx(i,j)**2+vely(i,j)**2)

			   if(i<=nx/2) then
	        
			      pre(i,j)=enin(i,j)*(gammal-1.)*den(i,j)
                    
	              sound(i,j)=sqrt(gammal*pre(i,j)/den(i,j))
			   else
           	      pre(i,j)=enin(i,j)*(gammar-1.)*den(i,j)
      
	              sound(i,j)=sqrt(gammar*pre(i,j)/den(i,j))
	           endif
          enddo
	enddo

	end

c------------------------------------------------------------------------------------
c                            getdt
c------------------------------------------------------------------------------------      
      subroutine getdt(sound,amass,tau,vertx,verty,Vstarx,Vstary,nx,ny,
     &	             dt0,dt)
      implicit double precision (a-h,o-z)

      double precision amass(1:nx,1:ny)
      double precision sound(1:nx,1:ny),tau(1:nx,1:ny)

	double precision vertx(0:nx,0:ny),verty(0:nx,0:ny)
      double precision Vstarx(0:nx,0:ny),Vstary(0:nx,0:ny)
     
	CE=0.3
	CV=0.1
	CM=1.01
	tE=9999.
	tV=9999.
      alanmda=9999.

	do i=1,nx
	   do j=1,ny
                            ! 寻找单元中最小的长度            
		  alanmda=999999.

	               xm=vertx(i-1,j-1)
	               ym=verty(i-1,j-1)
	               xn=vertx(i,j-1)
	               yn=verty(i,j-1)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离

	               xm=vertx(i-1,j-1)
	               ym=verty(i-1,j-1)
	               xn=vertx(i,j)
	               yn=verty(i,j)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离

	               xm=vertx(i-1,j-1)
	               ym=verty(i-1,j-1)
	               xn=vertx(i-1,j)
	               yn=verty(i-1,j)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离

	               xm=vertx(i,j-1)
	               ym=verty(i,j-1)
	               xn=vertx(i,j)
	               yn=verty(i,j)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离

	               xm=vertx(i,j-1)
	               ym=verty(i,j-1)
	               xn=vertx(i-1,j)
	               yn=verty(i-1,j)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离

	               xm=vertx(i,j)
	               ym=verty(i,j)
	               xn=vertx(i-1,j)
	               yn=verty(i-1,j)
	               call Length(xm,ym,xn,yn,alength)
	               alanmda=min(alanmda,alength)    !单元中边的最小距离        	
	         
             tE=min(CE*alanmda/sound(i,j),tE)
	   enddo            
	enddo 

      do i=1,nx
	   do j=1,ny

	      dAdt=0.
c1	          
			    ir=i-1
	            jr=j-1
	            ir1=i
	            jr1=j-1

	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            ur=Vstarx(ir,jr)    ! u star r
	            vr=Vstary(ir,jr)    ! v star r
	            ur1=Vstarx(ir1,jr1)
	            vr1=Vstary(ir1,jr1)
               
                dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
c2
			    ir=i
	            jr=j-1
	            ir1=i
	            jr1=j

	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            ur=Vstarx(ir,jr)    ! u star r
	            vr=Vstary(ir,jr)    ! v star r
	            ur1=Vstarx(ir1,jr1)
	            vr1=Vstary(ir1,jr1)
               
                dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
c3
			    ir=i
	            jr=j
	            ir1=i-1
	            jr1=j

	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            ur=Vstarx(ir,jr)    ! u star r
	            vr=Vstary(ir,jr)    ! v star r
	            ur1=Vstarx(ir1,jr1)
	            vr1=Vstary(ir1,jr1)
               
                dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
c4
			    ir=i-1
	            jr=j
	            ir1=i-1
	            jr1=j-1

	            xr=vertx(ir,jr)           
	            yr=verty(ir,jr)
	            xr1=vertx(ir1,jr1)
	            yr1=verty(ir1,jr1)

	            ur=Vstarx(ir,jr)    ! u star r
	            vr=Vstary(ir,jr)    ! v star r
	            ur1=Vstarx(ir1,jr1)
	            vr1=Vstary(ir1,jr1)
               
                dAdt=dAdt+(ur*yr1+vr1*xr-ur1*yr-vr*xr1)*0.5
	               
 		  	tV=min(CV*amass(i,j)*tau(i,j)/abs(dAdt),tV)   ! 这种方式可以考虑下  

	   enddo
	enddo
 
         dt=min(min(tV,tE),dt0*CM)

      end

c------------------------------------------------------------------------------------
c                            Length(r,r+1)
c------------------------------------------------------------------------------------
      subroutine Length(x1,y1,x2,y2,alength)
      implicit double precision (a-h,o-z)
      
      tt=(x1-x2)**2+(y1-y2)**2
	if(tt<0.) then
	pause
	endif
      alength=sqrt((x1-x2)**2+(y1-y2)**2)      

      end 
c------------------------------------------------------------------------------------
c                            Normal  N^{k-1}_{k}
c------------------------------------------------------------------------------------
      subroutine Normal(x1,y1,x2,y2,anormalx,anormaly)
      implicit double precision (a-h,o-z)
      
      tt=(x1-x2)**2+(y1-y2)**2

      if(tt<0.) then
	pause
	endif

	tempt=sqrt((x1-x2)**2+(y1-y2)**2)

      anormalx=(y1-y2)/tempt
	anormaly=(x2-x1)/tempt

      end

c------------------------------------------------------------------------------------
c                 three   求解 3×3 行列式
c------------------------------------------------------------------------------------
      subroutine three(a11,a12,a13,a21,a22,a23,a31,a32,a33,det)
      implicit double precision (a-h,o-z)

      det=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)
     &   +a13*(a21*a32-a22*a31)
      
	end
	