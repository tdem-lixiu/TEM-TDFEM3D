    SUBROUTINE GRID(XYZ,I3,N3,NI,NS,WE,SX)
    IMPLICIT NONE
    INTEGER I,J,K,NS,WE,SX,N,M,L,K1,KK
    !N为计算区域的单元数，M为棱边数，L为节点数
    INTEGER I3((NS-1)*(WE-1)*(SX-1),8),N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2)
    REAL(8) XYZ(NS*WE*SX,3),SUM,DX(WE),DY(NS),DZ(SX)
    N=(NS-1)*(WE-1)*(SX-1)
    L=NS*WE*SX
    M=(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
    !得到节点编号
      !DIAN
      I3(1,1)=1
      I3(1,2)=2
      I3(1,3)=WE+2
      I3(1,4)=WE+1      
      I3(1,5)=I3(1,1)+WE*NS
      I3(1,6)=I3(1,2)+WE*NS
      I3(1,7)=I3(1,3)+WE*NS
      I3(1,8)=I3(1,4)+WE*NS
      !XIAN
      DO I=2,WE-1
      DO J=1,8
      I3(I,J)=I3(I-1,J)+1
      ENDDO
      ENDDO
      !MIAN
      K=WE
      DO I=2,NS-1
      DO K1=1,WE-1
      DO J=1,8
      I3(K,J)=I3(K-WE+1,J)+WE
      ENDDO
      K=K+1
      ENDDO
      ENDDO
      !TI
      DO I=2,(SX-1)
      DO K1=1,(WE-1)*(NS-1)
      DO J=1,8
      I3(K,J)=I3(K-(WE-1)*(NS-1),J)+WE*NS
      ENDDO
      K=K+1
      ENDDO
      ENDDO 
      !得到棱边编号
      !DIAN
      N3(1,1)=1
      N3(1,2)=WE-1+1
      N3(1,3)=N3(1,1)+(WE-1)*NS
      N3(1,4)=N3(1,2)+(WE-1)*NS
      N3(1,5)=1+(WE-1)*SX*NS
      N3(1,6)=N3(1,5)+(NS-1)*WE
      N3(1,7)=N3(1,5)+1
      N3(1,8)=N3(1,7)+(NS-1)*WE
      N3(1,9)=1+(WE-1)*NS*SX+(NS-1)*WE*SX
      N3(1,10)=N3(1,9)+1
      N3(1,11)=N3(1,9)+WE
      N3(1,12)=N3(1,11)+1
      !XIAN
      DO I=2,WE-1
      DO J=1,12
      N3(I,J)=N3(I-1,J)+1
      ENDDO
      ENDDO
      !MIAN
      K=WE
      DO I=2,NS-1
      DO K1=1,WE-1
      DO J=1,4
      N3(K,J)=N3(K-WE+1,J)+WE-1
      ENDDO
      DO J=5,12
      N3(K,J)=N3(K-WE+1,J)+WE
      ENDDO
      K=K+1
      ENDDO
      ENDDO
      !TI
      DO I=2,(SX-1)
      DO K1=1,(WE-1)*(NS-1)
      DO J=1,4
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+(WE-1)*NS
      ENDDO
      DO J=5,8
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+WE*(NS-1)
      ENDDO
      DO J=9,12
      N3(K,J)=N3(K-(WE-1)*(NS-1),J)+WE*NS
      ENDDO
      K=K+1
      ENDDO
      ENDDO 
      !得到棱边对应的节点
      !X
      DO I=1,WE-1
      NI(I,1)=I
      NI(I,2)=I+1
      ENDDO
      DO I=WE,(WE-1)*NS
      NI(I,1)=NI(I-(WE-1),1)+WE
      NI(I,2)=NI(I-(WE-1),2)+WE
      ENDDO
      DO I=(WE-1)*NS+1,(WE-1)*NS*SX
      NI(I,1)=NI(I-(WE-1)*NS,1)+WE*NS
      NI(I,2)=NI(I-(WE-1)*NS,2)+WE*NS
      ENDDO 
      !Y
      K1=1
      DO I=(WE-1)*NS*SX+1,(WE-1)*NS*SX+WE
      NI(I,1)=K1
      NI(I,2)=K1+WE
      K1=K1+1
      ENDDO
      DO I=(WE-1)*NS*SX+WE+1,(WE-1)*NS*SX+(NS-1)*WE
      NI(I,1)=NI(I-WE,1)+WE
      NI(I,2)=NI(I-WE,2)+WE
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE+1,(WE-1)*NS*SX+(NS-1)*WE*SX
      NI(I,1)=NI(I-WE*(NS-1),1)+WE*NS
      NI(I,2)=NI(I-WE*(NS-1),2)+WE*NS
      ENDDO
      !Z
      K1=1
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE
      NI(I,1)=K1
      NI(I,2)=K1+WE*NS
      K1=K1+1
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+WE+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS
      NI(I,1)=NI(I-WE,1)+WE
      NI(I,2)=NI(I-WE,2)+WE
      ENDDO
      DO I=(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS+1,(WE-1)*NS*SX+(NS-1)*WE*SX+WE*NS*(SX-1)
      NI(I,1)=NI(I-WE*NS,1)+WE*NS
      NI(I,2)=NI(I-WE*NS,2)+WE*NS
      ENDDO
      
      
      
      
      
      !得到节点坐标
       OPEN(100,FILE='xyz.txt')
       READ(100,*)
       READ(100,*)
       READ(100,*)(DX(I),I=1,WE)    !X方向节点坐标
       READ(100,*)
       READ(100,*)(DY(I),I=1,NS)    !Y方向节点坐标
       READ(100,*)
       READ(100,*)(DZ(I),I=1,SX)    !Z方向节点坐标
       CLOSE(100)
       DO I=1,WE
       DO J=1,NS
       DO K=1,SX
       XYZ((K-1)*WE*NS+(J-1)*WE+I,1)=DX(I)
       XYZ((K-1)*WE*NS+(J-1)*WE+I,2)=DY(J)
       XYZ((K-1)*WE*NS+(J-1)*WE+I,3)=DZ(K)
       ENDDO
       ENDDO
       ENDDO
      
    ENDSUBROUTINE