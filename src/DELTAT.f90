
    
    !
    
    
    
    !Դ�ļ��أ��������ص����β�
    SUBROUTINE YUAN11(L,I0,PI,WE,NS,SX,M,J)
    IMPLICIT NONE
    INTEGER M,WE,NS,SX,I,N,L
    REAL(8) T,PI,U,J((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),I0,DELTAT
    IF (M==1) THEN
        U=I0/1.0E-6/100.0    !������
    ENDIF
    IF (M>=2.AND.M<=24) THEN
        U=0                   !����ʱ��
    ENDIF
    IF (M==25) THEN
        U=-I0/1.0E-6/100.0    !�½���
    ENDIF
    IF (M>=26.AND.M<=39) THEN
        U=0                   !�ضϺ�
    ENDIF 
     J=0
     
     !Դ��λ�ü�����
    DO I=22*(WE-1)*NS+((NS+1)/2.0-L/10.0/2.0)*(WE-1)-(WE-1)/2.-L/10.0/2.0+1.0,22*(WE-1)*NS+((NS+1)/2.0-L/10.0/2.0)*(WE-1)-(WE-1)/2.+L/10.0/2.0   !  Դ��λ��
        J(I)=-U
        J(I+L/10.0*(WE-1))=U
    ENDDO
    DO I=1,L/10.0
        J((WE-1)*NS*SX+22*WE*(NS-1)+((NS-1)/2.-L/10./2.)*WE+(WE+1)/2.-L/10./2+(I-1)*WE)=U
        J((WE-1)*NS*SX+22*WE*(NS-1)+((NS-1)/2.-L/10./2.)*WE+(WE+1)/2.-L/10./2+(I-1)*WE+L/10.0)=-U
    ENDDO
    ENDSUBROUTINE
    
    