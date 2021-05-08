
    SUBROUTINE CALCULATE11(XYZ,NI,N3,I3,PP,NS,WE,SX,N)
    IMPLICIT NONE
    INTEGER I,N,WE,NS,SX,KK,III,II,JJ,k,L
    INTEGER I3((NS-1)*(WE-1)*(SX-1),8),N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2)
    REAL(8) E1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),DELTAT(N),PI,J((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    REAL(8) T(N),U0,PP((NS-1)*(WE-1)*(SX-1)),XYZ(WE*NS*SX,3)
    REAL(8) AA((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),I0,EE((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    REAL(8) B((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),E2((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    REAL(8) HHH(WE*NS*SX,3),HH2(WE*NS*SX,3),X((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    INTEGER XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
    PARAMETER(XMIN=28,XMAX=36,YMIN=28,YMAX=36,ZMIN=23,ZMAX=23)   !输出的范围
    INTEGER RecFile1(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1),RecFile(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1)
    CHARACTER(80) FILENAME(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1),FILENAME1(XMAX-XMIN+1,YMAX-YMIN+1,ZMAX-ZMIN+1)
    CHARACTER(2)STRX,STRY,STRZ
    I0=100.0            !电流
    PI=3.14159265358979323846
    U0=4*PI*1.0E-7
    L=100                                        !L为回线源边长
    OPEN(130,FILE='TT.txt')
    DO II=1,N
       READ(130,*)III,T(II)          !读取时间道
    ENDDO
    DO II=1,n/100.
        READ(130,*)III,DELTAT(II)    !读取时间间隔
    ENDDO
    CLOSE(130)
    
     !--------输出-------
    DO KK=XMIN,XMAX,1
        DO JJ=YMIN,YMAX,1
            DO K=ZMIN,ZMAX,1
                RecFile(KK,JJ,K)=20000+200*KK+JJ
                RecFile1(KK,JJ,K)=20000+200*KK+JJ+2000000
            ENDDO
        ENDDO
    ENDDO
    DO KK=XMIN,XMAX,1
        DO JJ=YMIN,YMAX,1
            DO K=ZMIN,ZMAX,1
            WRITE(STRX,'(I2.2)')KK
            STRX=TRIM(ADJUSTl(STRX))
            STRX=TRIM(ADJUSTR(STRX))
            WRITE(STRY,'(I2.2)')JJ
            STRY=TRIM(ADJUSTl(STRY))
            STRY=TRIM(ADJUSTR(STRY))
            WRITE(STRZ,'(I2.2)')K
            STRZ=TRIM(ADJUSTl(STRZ))
            STRZ=TRIM(ADJUSTR(STRZ))
            FILENAME(KK,JJ,K)='E'//STRX//'-'//STRY//'-'//STRZ//'.DAT' 
            FILENAME1(KK,JJ,K)='Decay-Curve-'//STRX//'-'//STRY//'-'//STRZ//'.DAT' 
            OPEN(RecFile(KK,JJ,K),FILE=FILENAME(KK,JJ,K)) 
            OPEN(RecFile1(KK,JJ,K),FILE=FILENAME1(KK,JJ,K))
            ENDDO
        ENDDO
    ENDDO
    
    !--------------
    E2=0
    DO I=1,N/100.
        PRINT*,I
        
    CALL HECHENG(FILENAME,FILENAME1,RecFile,RecFile1,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,E2,I,I3,N3,XYZ,PP,U0,WE,NS,SX,DELTAT(I),T,L,N,PI,I0,NI)
    
     
    ENDDO
    DO KK=XMIN,XMAX,1
        DO JJ=YMIN,YMAX,1
            DO K=ZMIN,ZMAX,1
              CLOSE(RecFile(KK,JJ,K))
              CLOSE(RecFile1(KK,JJ,K))
            ENDDO
        ENDDO
     ENDDO
    ENDSUBROUTINE
    
    
    
    

    
   