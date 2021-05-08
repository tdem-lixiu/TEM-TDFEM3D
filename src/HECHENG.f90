     
    !该子程序为单元合成
    SUBROUTINE HECHENG(FILENAME,FILENAME1,RecFile,RecFile1,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,E2,I1,I3,N3,XYZ,PP,U0,WE,NS,SX,DELTAT,T,LL,NN,PI,I0,NI)
                                                                
      IMPLICIT NONE
      INTEGER NS,WE,SX,L,M,I,J,K1,K2,MM
      
      INTEGER XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,I1,NN,LL,n,K,KK,JJ,ND
      REAL(8) T(NN),PI,I0,J1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      INTEGER RecFile1(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1),RecFile(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1)
      CHARACTER(80) FILENAME(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1),FILENAME1(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1)
      REAL(8) HHH(WE*NS*SX,3),HH2(WE*NS*SX,3)   !HHH为每个节点的电场，HH2为每个节点的衰减电压
     
     
      
      INTEGER N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2),I3((NS-1)*(WE-1)*(SX-1),8)
      INTEGER IA((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)+1),N1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      REAL(8) XYZ(WE*NS*SX,3),PP((NS-1)*(WE-1)*(SX-1)),DELTAT,B((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      REAL(8) X((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),E2((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      REAL(8) EA(12,12),U0,E1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
      REAL(8),ALLOCATABLE:: AA(:)
      INTEGER,ALLOCATABLE:: JA(:)
      
      TYPE NODE
        REAL(8) A
        INTEGER KI,KJ
        TYPE(NODE),POINTER:: NEXT
      END TYPE NODE
      
      TYPE NODE1
        TYPE(NODE),POINTER:: HEAD
      END TYPE NODE1   
         
      TYPE(NODE1),ALLOCATABLE:: H(:)
      TYPE(NODE),POINTER:: Q,P,P1,Q1
      
      
      INTEGER*8 pt(64)
     INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl 
     integer,allocatable::perm(:)
     INTEGER iparm(64)
     Integer mkl_get_max_threads
     external mkl_get_max_threads
     allocate(perm((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)))
     
     ALLOCATE(H((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)))
     
      
      DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
       NULLIFY(H(I).HEAD)
      ENDDO
      M=0
      DO L=1,(NS-1)*(WE-1)*(SX-1)
         
       CALL UKE1(DELTAT,L,XYZ,N3,NI,I3,NS,SX,WE,EA,PP,U0)

        DO I=1,12
        DO J=I,12
            IF(ABS(EA(I,J))>0.0) THEN
                IF(N3(L,I)<=N3(L,J))THEN
                K1=N3(L,I)
                K2=N3(L,J)
                ELSE
                K1=N3(L,J)
                K2=N3(L,I)
                ENDIF
                
                
       ! 判断是否为空表头
                IF(.NOT.ASSOCIATED(H(K1).HEAD)) THEN
                    ALLOCATE(Q)
                    Q.A=EA(I,J)
                    Q.KI=K1
                    Q.KJ=K2
                    NULLIFY(Q.NEXT) 
                  H(K1).HEAD=>Q 
                  M=M+1
                  N1(K1)=1
                ELSE
                    P1=>H(K1).HEAD
                    NULLIFY(P)
                  !寻找数组中元素的相同位置
                 WW: DO WHILE(ASSOCIATED(P1))
                       IF(K2==P1.KJ)THEN
                         P=>P1
                        EXIT WW
                       ELSE
                         P1=>P1.NEXT     !链表的后退方向
                       ENDIF
                     ENDDO WW
                  !如果有重复，进行累加，如果没有存到后一个节点
                     IF(ASSOCIATED(P))THEN
                       P.A=P.A+EA(I,J)
                     ELSE
                          ALLOCATE(Q)
                           Q.A=EA(I,J)
                           Q.KI=K1
                           Q.KJ=K2
                           NULLIFY(Q.NEXT)
                         IF(Q.KJ<=H(K1).HEAD.KJ)THEN
                       Q.NEXT=>H(K1).HEAD
                       H(K1).HEAD=>Q
                       M=M+1
                     N1(K1)=N1(K1)+1    
                     ELSE
                       NULLIFY(P1) 
                       P1=>H(K1).HEAD
                  WW2:     DO WHILE(ASSOCIATED(P1).AND.Q.KJ>P1.KJ)
                             Q1=>P1
                            P1=>P1.NEXT  
                            EXIT WW2             
                          ENDDO WW2
                       IF(ASSOCIATED(P1))THEN 
                       M=M+1
                       N1(K1)=N1(K1)+1
                        Q.NEXT=>Q1.NEXT
                        Q1.NEXT=>Q
                       ELSE
                       M=M+1
                       N1(K1)=N1(K1)+1
                       Q1.NEXT=>Q
                       ENDIF
                     ENDIF
                     ENDIF
               ENDIF  
            ENDIF 
           ENDDO
          ENDDO
      ENDDO  
    ! 转化为静态数组  !释放链表
       ALLOCATE(JA(M),AA(M))
       IA(1)=1
       J=1
       DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
          P=>H(I).HEAD
          DO WHILE(ASSOCIATED(P))
             AA(J)=P.A
             JA(J)=P.KJ
             Q1=>P.NEXT
             DEALLOCATE(P)
             P=>Q1
             J=J+1 
          ENDDO
             IA(I+1)=J
             IF(N1(I)/=(IA(I+1)-IA(I)).OR.N1(I)==0) THEN
             PRINT*,'YOUCUOWU'
             PRINT*,N1(I),IA(I+1)-IA(I),I
             PAUSE
             ENDIF
       ENDDO
       DEALLOCATE(H)
       ND=(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
       
       
    !****************解方程***************************
    
    !调用pardiso求解器求解方程组，详情请参考MKL函数库使用说明
    
    pt=0 !指针初始化
    maxfct=1;mnum=1;mtype=2;n=ND;perm=0;nrhs=1;x=0.0
    iparm(1) = 0 ! iparm使用默认值
    iparm(3) = mkl_get_max_threads()  ! iparm（3）并行线程数由mkl获得，一般有几个cpu就有几个线程
    error = 0 ! 错误信息初值?
    msglvl=0
        
        E1=E2
        CALL YUAN11(LL,I0,PI,WE,NS,SX,I1,J1)   !源的加载
        
        CALL YOUHECHENG(U0,PP,WE,NS,SX,XYZ,I3,N3,E1,J1,DELTAT,B)   !右端项的单元合成
       
        !**********解方程**************
    
    
    phase=12
    call pardiso (pt, maxfct, mnum, mtype, phase, n, aa, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    print*,'fenjie'
    if(error/=0) THEN
    print*,'SHIBAI',error
    PAUSE
    ENDIF
    
    
   
    phase=33
    call pardiso (pt, maxfct, mnum, mtype, phase, n, aa, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    print*,'jiefang'
    if(error/=0) THEN
    print*,'shibai1',error
    PAUSE
    ENDIF
    E2=X
    
    
    !****************************方程求解结束************************
    
    
    
    do MM=(I1-1)*100+2,(I1-1)*100+100
        PRINT*,'zong=',I1,MM
        E1=E2
        
        CALL YOUHECHENG(U0,PP,WE,NS,SX,XYZ,I3,N3,E1,J1,DELTAT,B)   !右端项的合成
        PRINT*,'S'
        !************************解方程********************
        phase=33
    call pardiso (pt, maxfct, mnum, mtype, phase, n, aa, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    print*,'jiefang'
    if(error/=0) THEN
    print*,'shibai3',error
    PAUSE
    ENDIF
    E2=X
    if (MM>2700.and.modulo(MM,3)==0) then
    CALL NB(E2,HHH,XYZ,WE,NS,SX,NI)    !棱边转换
    CALL HZA(HHH,HH2,XYZ,WE,NS,SX)     !电场转衰减电压
    DO KK=XMIN,XMAX,1
        DO JJ=YMIN,YMAX,1
            DO K=ZMIN,ZMAX,1
                !电场
              WRITE(RecFile(KK,JJ,K),'(7E)')XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,1),XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,2),XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,3),T(MM),HHH(KK+(JJ-1)*WE+(K-1)*WE*NS,1),HHH(KK+(JJ-1)*WE+(K-1)*WE*NS,2),HHH(KK+(JJ-1)*WE+(K-1)*WE*NS,3)
                !衰减电压
              WRITE(RecFile1(KK,JJ,K),'(7E)')XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,1),XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,2),XYZ(KK+(JJ-1)*WE+(K-1)*WE*NS,3),T(MM),HH2(KK+(JJ-1)*WE+(K-1)*WE*NS,1),HH2(KK+(JJ-1)*WE+(K-1)*WE*NS,2),HH2(KK+(JJ-1)*WE+(K-1)*WE*NS,3)
           ENDDO
        ENDDO
     ENDDO
    endif
    enddo
    !**************************方程求解结束*****************
   
    phase=-1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, aa, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    deallocate(perm)
    DEALLOCATE(JA,AA)
    !************************方程求解结束**************
    
       !PRINT*,'FANGCHENGQIUJIEZHONG'
       !PRINT*,DD
       !CALL PARDISO0(DD,AA,JA,M,B,IA,ND,E2)
      
      
    ENDSUBROUTINE
      !***************************************
    
      !该子程序为单元分析
      SUBROUTINE UKE1(DELTAT,K,XYZ,N3,NI,I3,NS,SX,WE,EA,PP,U0)
      IMPLICIT NONE
      INTEGER I,J,K,NS,SX,WE
      INTEGER N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2),I3((NS-1)*(WE-1)*(SX-1),8)
      REAL(8) A1(4,4),A2(4,4),A3(4,4),A4(4,4),PP((NS-1)*(WE-1)*(SX-1)),A31(4,4),LX,LY,LZ,U0,XYZ(NS*WE*SX,3)
      REAL(8) EA(12,12),EA1(12,12),EA2(12,12),DELTAT
      DATA A1/2.,-2.,1.,-1.,-2.,2.,-1.,1.,1.,-1.,2.,-2.,-1.,1.,-2.,2./
      DATA A2/2.,1.,-2.,-1.,1.,2.,-1.,-2.,-2.,-1.,2.,1.,-1.,-2.,1.,2./
      DATA A3/2.,-2.,1.,-1.,1.,-1.,2.,-2.,-2.,2.,-1.,1.,-1.,1.,-2.,2./
      DATA A4/4.,2.,2.,1.,2.,4.,1.,2.,2.,1.,4.,2.,1.,2.,2.,4./
      
      
      A31=TRANSPOSE(A3)
      LX=ABS(XYZ(I3(K,1),1)-XYZ(I3(K,2),1))
      LY=ABS(XYZ(I3(K,1),2)-XYZ(I3(K,4),2))
      LZ=ABS(XYZ(I3(K,1),3)-XYZ(I3(K,5),3))
      DO I=1,4
      DO J=1,4
      EA1(I,J)=LX*LZ/6./LY*A1(I,J)+LX*LY/6./LZ*A2(I,J)
      EA1(I+4,J+4)=LX*LY/6./LZ*A1(I,J)+LZ*LY/6./LX*A2(I,J)
      EA1(I+8,J+8)=LY*LZ/6./LX*A1(I,J)+LX*LZ/6./LY*A2(I,J)
      EA1(I+4,J)=-LZ/6.*A31(I,J)
      EA1(I+8,J)=-LY/6.*A3(I,J)
      EA1(I+8,J+4)=-LX/6.*A31(I,J)
      
      EA1(I,J+4)=-LZ/6.*A3(I,J)   
      EA1(I,J+8)=-LY/6.*A31(I,J)
      EA1(I+4,J+8)=-LX/6.*A3(I,J)
      
      ENDDO
      ENDDO
      EA2=0.0
      DO I=1,4
      DO J=1,4
      EA2(I,J)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)/DELTAT                        !!上给出DELTAT
      EA2(I+4,J+4)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)/DELTAT
      EA2(I+8,J+8)=U0*PP(K)*LX*LY*LZ/36.0*A4(I,J)/DELTAT
      ENDDO
      ENDDO 
      DO I=1,12
      DO J=1,12
      EA(I,J)=EA1(I,J)+EA2(I,J)
      ENDDO
      ENDDO
      ENDSUBROUTINE
