     !右端项的合成
    SUBROUTINE YOUHECHENG(U0,PP,WE,NS,SX,XYZ,I3,N3,E1,J,DELTAT,B)
    
    IMPLICIT NONE
    INTEGER L,K1,WE,NS,SX
    INTEGER N3((NS-1)*(WE-1)*(SX-1),12),I3((NS-1)*(WE-1)*(SX-1),8)
    REAL(8) U0,PP((NS-1)*(WE-1)*(SX-1)),XYZ(WE*NS*SX,3)
    REAL(8) E1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),J((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    REAL(8) DELTAT,B((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),YD(12)
    
    B=0
    DO L=1,(WE-1)*(NS-1)*(SX-1)
        CALL FENXI(U0,PP,WE,NS,SX,XYZ,I3,N3,L,DELTAT,E1,J,YD)   !单元分析
    DO K1=1,12
    B(N3(L,K1))=B(N3(L,K1))+YD(K1)
    ENDDO
    
    ENDDO
    ENDSUBROUTINE
    
    
    
    
    
    !!右端项单元分析
    SUBROUTINE FENXI(U0,PP,WE,NS,SX,XYZ,I3,N3,I,DELTAT,E1,J,YD)                     !!DELTAT 和 J 传进来
    IMPLICIT NONE
    INTEGER I,WE,NS,SX,K1,K2,K
    INTEGER N3((NS-1)*(WE-1)*(SX-1),12),I3((NS-1)*(WE-1)*(SX-1),8)
    REAL(8) U0,PP((NS-1)*(WE-1)*(SX-1)),XYZ(WE*NS*SX,3),EZ(12,12)
    REAL(8) E1((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),J((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1))
    REAL(8) DELTAT,YD(12),LX,LY,LZ,A4(4,4),EB(12)
    REAL(8) YB(12)
    
    DATA A4/4.,2.,2.,1.,2.,4.,1.,2.,2.,1.,4.,2.,1.,2.,2.,4./
    !DO I=1,(WE-1)*(NS-1)*(SX-1)
    LX=ABS(XYZ(I3(I,1),1)-XYZ(I3(I,2),1))
    LY=ABS(XYZ(I3(I,1),2)-XYZ(I3(I,4),2))
    LZ=ABS(XYZ(I3(I,1),3)-XYZ(I3(I,5),3))
    
    EZ=0.0
    DO K1=1,4
    DO K2=1,4
    EZ(K1,K2)=U0*PP(I)*LX*LY*LZ/36.0*A4(K1,K2)/DELTAT
    EZ(K1+4,K2+4)=U0*PP(I)*LX*LY*LZ/36.0*A4(K1,K2)/DELTAT
    EZ(K1+8,K2+8)=U0*PP(I)*LX*LY*LZ/36.0*A4(K1,K2)/DELTAT
    ENDDO
    ENDDO
    
    DO K1=1,12
       EB(K1)=E1(N3(I,K1))
    ENDDO   
    EB=MATMUL(EZ,EB)                   !!电场单元分析
    
    !!源的单元分析----------------------------
    
       EZ=0.0
    DO K1=1,4
    DO K2=1,4
    EZ(K1,K2)=U0*LX*LY*LZ/36.0*A4(K1,K2)
    EZ(K1+4,K2+4)=U0*LX*LY*LZ/36.0*A4(K1,K2)
    EZ(K1+8,K2+8)=U0*LX*LY*LZ/36.0*A4(K1,K2)
    ENDDO
    ENDDO
    
    
    DO K1=1,12
       YB(K1)=J(N3(I,K1))
    ENDDO   
    YB=MATMUL(EZ,YB)                      !!源的单元分析
    
    !!右端两项之和---------------
    
    YD=EB-YB                             !!右端两项只和
   
   !!单元合成 
    
    !DO K1=1,12
    !B(N3(I,K1))=B(N3(I,K1))+U0*EB(K1)
    !ENDDO
    
    
    
    
    !ENDDO
    ENDSUBROUTINE
    