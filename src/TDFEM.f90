
    
    !�ó���Ϊ��ʱ����ʸ������Ԫ˲������ά���ݡ��ļ������   

 !�˳���Ľṹ����Ҫ�ӳ���Ĺ������£�
 
 !�˳�����Ҫ��Ϊ���󲿷֣��ֱ�Ϊ���ݶ�ȡ�����ݼ�������������

 !���������ݶ�ȡ��������Ҫ���������ӳ���,
 !��һ�� GET_PP(PP,NS,SX,WE)�� �˳���Ĺ����Ƕ�ȡÿ����Ԫ�ĵ����ʣ�����PP��ʾ
 !ÿ����Ԫ�ĵ����ʣ�WE,NS,SX��ʾģ����X,Y,Z���������ϵĽڵ������ڶ��� 
 !GRID(XYZ,I3,N3,NI,NS,WE,SX)���˳���Ĺ����Ƕ�ģ�͵�ÿ����Ԫ����߼��ڵ����
 !��ţ����Ҷ�ȡÿ���ڵ�����꣬����XYZ����ÿ���ڵ�����꣬I3����ÿ����Ԫ��8
 !���ڵ㣬N3����ÿ����Ԫ��12����ߣ�NI����ÿ����ߵ�2���ڵ㡣

 !����ڼ����������Ҫ�����ĸ��ӳ���
 !��һ�� UKE1(DELTAT,K,XYZ,N3,NI,I3,NS,SX,WE,EA,PP,U0)���˳���Ĺ����ǽ��е�Ԫ
 !����������DELTAT��ʾʱ�䲽����K��ʾ��Ԫ�ı�ţ�EA��ʾ��Ԫ����U0��ʾ�ŵ�
 !�ʣ��ڶ��� HECHENG()���˳���Ĺ����ǵ�Ԫ�ϳɣ������� YUAN11()���˳���Ĺ���
 !��Դ�ļ��أ������� YOUHECHENG()���˳���Ĺ�������ⷽ������Ҷ��

 !��������������������Ҫ���������ӳ���
 !��һ�� NB(E2,HHH,XYZ,WE,NS,SX,NI)���˳���Ĺ��������ת�����������ɶ������
 !���ĵ�ĵ糡ֵת�����ڵ��ϣ�����E2��ʾ��������ϵĵ糡ֵ��HHH��ʾת����ĸ�
 !���ڵ��ϵĵ糡ֵ���ڶ��� HZA(HHH,HH2,XYZ,WE,NS,SX)���˳���Ĺ����ǽ��ڵ��ϵ�
 !�糡ת��ΪdB/dt������HH2��ʾת����ĸ����ڵ��ϵ�dB/dt��
    
 PROGRAM FEM  
 IMPLICIT NONE
 INTEGER WE,NS,SX,I,J,N,nd                  
 PARAMETER(WE=63,NS=63,SX=53,N=3900)!WE��NS,SX, X,Y,Z����ڵ�����N��ʾʱ�����
 INTEGER I3((NS-1)*(WE-1)*(SX-1),8),N3((NS-1)*(WE-1)*(SX-1),12),NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2)
 !!I3ÿ����Ԫ��8����  N3ÿ����Ԫ��12����  NIÿ����ߵ�2���ڵ�
 REAL(8) XYZ(NS*WE*SX,3),PP((NS-1)*(WE-1)*(SX-1))  !!XYZÿ���ڵ�����꣬PPÿ����Ԫ�ĵ�����
 
 CALL GET_PP(PP,NS,SX,WE)!���ӳ����ȡ��ÿ����Ԫ�ĵ�����
 CALL GRID(XYZ,I3,N3,NI,NS,WE,SX)!���ӳ����ÿ����Ԫ����߼��ڵ�����˱��
 
 CALL CALCULATE11(XYZ,NI,N3,I3,PP,NS,WE,SX,N)!���ӳ���Ϊ�������
    ENDPROGRAM
    
   
    !���ӳ���Ϊ���ת��
    SUBROUTINE NB(HXYZ,HHH,XYZ,WE,NS,SX,NI)
    IMPLICIT NONE
    INTEGER WE,NS,SX,N,I,J,K,K1
    INTEGER NI((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),2)
    REAL(8) XYZ(NS*WE*SX,3),XY((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1),3)
    REAL(8) HXYZ((WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)),HHH(WE*NS*SX,3)
    REAL(8) X(WE-1),U(WE-1),X1(WE-2),U1(WE-2)
    REAL(8) XX(NS-1),UU(NS-1),XX1(NS-2),UU1(NS-2)
    REAL(8) XXX(SX-1),UUU(SX-1),XXX1(SX-2),UUU1(SX-2)
    HHH=0.0
    DO I=1,(WE-1)*NS*SX+WE*(NS-1)*SX+WE*NS*(SX-1)
    XY(I,1)=(XYZ(NI(I,1),1)+XYZ(NI(I,2),1))/2.0
    XY(I,2)=(XYZ(NI(I,1),2)+XYZ(NI(I,2),2))/2.0
    XY(I,3)=(XYZ(NI(I,1),3)+XYZ(NI(I,2),3))/2.0
    ENDDO
      
    !EX
    DO I=1,SX*NS
        DO J=1,WE-1
           X(J)=XY(J+(I-1)*(WE-1),1)
           U(J)=HXYZ(J+(I-1)*(WE-1))
        ENDDO
        DO J=1,WE-2
           X1(J)=XYZ(J+(I-1)*WE+1,1)
        ENDDO
        CALL SPLN3(X,U,WE-1,X1,U1,WE-2)
        DO J=1,WE-2
           HHH(J+(I-1)*WE+1,1)=U1(J)
        ENDDO
    ENDDO
     !EY
    DO I=1,SX
        DO K=1,WE
            DO J=1,NS-1
               XX(J)=XY((WE-1)*NS*SX+(I-1)*WE*(NS-1)+K+(J-1)*WE,2)
               UU(J)=HXYZ((WE-1)*NS*SX+(I-1)*WE*(NS-1)+K+(J-1)*WE)
            ENDDO
            DO J=1,NS-2
               XX1(J)=XYZ((I-1)*WE*NS+K+(J-1)*WE+WE,2)
            ENDDO
            CALL SPLN3(XX,UU,NS-1,XX1,UU1,NS-2)
            DO J=1,NS-2
               HHH((I-1)*WE*NS+K+(J-1)*WE+WE,2)=UU1(J)
            ENDDO
        ENDDO
    ENDDO
     !EZ
    DO I=1,WE*NS
        DO J=1,SX-1
           XXX(J)=XY((WE-1)*NS*SX+WE*(NS-1)*SX+I+(J-1)*WE*NS,3)
           UUU(J)=HXYZ((WE-1)*NS*SX+WE*(NS-1)*SX+I+(J-1)*WE*NS)
        ENDDO
        DO J=1,SX-2
           XXX1(J)=XYZ(I+J*WE*NS,3)
        ENDDO
        CALL SPLN3(XXX,UUU,SX-1,XXX1,UUU1,SX-2)
        DO J=1,SX-2
           HHH(I+J*WE*NS,3)=UUU1(J)
        ENDDO
    ENDDO
    
    ENDSUBROUTINE
    !*********************************************
     SUBROUTINE SPLN3(X,Y,N,XX,YY,M)
      IMPLICIT NONE
      INTEGER J,N,M,I
      REAL(8) X(N),Y(N),DY(N),H(N),XX(M),YY(M)
      REAL(8) H0,H1,BETA,ALPHA,DY1,DYN,t1,t2,t3,t4,t5,t6,a,b
      t1=x(2)-x(1);t2=x(3)-x(1);t3=y(2)-y(1);t4=y(3)-y(1)
      t5=x(3)*x(3)-x(1)*x(1);t6=x(2)*x(2)-x(1)*x(1)
      a=(t4*t1-t3*t2)/(t1*t5-t2*t6);b=(t4*t6-t3*t5)/(t2*t6-t5*t1)
      DY1=2.0*a*x(1)+b
      t1=x(n-1)-x(n-2);t2=x(n)-x(n-2);t3=y(n-1)-y(n-2);t4=y(n)-y(n-2)
      t5=x(n)*x(n)-x(n-2)*x(n-2);t6=x(n-1)*x(n-1)-x(n-2)*x(n-2)
      a=(t4*t1-t3*t2)/(t1*t5-t2*t6);b=(t4*t6-t3*t5)/(t2*t6-t5*t1)
      DYN=2.0*a*x(n)+b
      DY(1)=0.0
      H(1)=DY1
      H0=X(2)-X(1)
      DO J=2,N-1
      H1=X(J+1)-X(J)
      ALPHA=H0/(H0+H1)
      BETA=(1.0-ALPHA)*(Y(J)-Y(J-1))/H0
      BETA=3.0*(BETA+ALPHA*(Y(J+1)-Y(J))/H1)
      DY(J)=-ALPHA/(2.0+(1.0-ALPHA)*DY(J-1))
      H(J)=(BETA-(1.0-ALPHA)*H(J-1))
      H(J)=H(J)/(2.0+(1.0-ALPHA)*DY(J-1))
      H0=H1
      ENDDO
      DY(N)=DYN
      DO J=N-1,1,-1
      DY(J)=DY(J)*DY(J+1)+H(J)
      ENDDO
      DO J=1,N-1
      H(J)=X(J+1)-X(J)
      ENDDO
      DO J=1,M
      IF (XX(J).GE.X(N)) THEN
      I=N-1
      ELSE
      I=1
60    IF (XX(J).GT.X(I+1)) THEN
      I=I+1
      GOTO 60
      END IF
      END IF
      H1=(X(I+1)-XX(J))/H(I)
      YY(J)=(3.0*H1*H1-2.0*H1*H1*H1)*Y(I)
      YY(J)=YY(J)+H(I)*(H1*H1-H1*H1*H1)*DY(I)
      H1=(XX(J)-X(I))/H(I)
      YY(J)=YY(J)+(3.0*H1*H1-2.0*H1*H1*H1)*Y(I+1)
      YY(J)=YY(J)-H(I)*(H1*H1-H1*H1*H1)*DY(I+1)
      ENDDO
      RETURN
    ENDSUBROUTINE
    !***********************************************
    !���ӳ��򽫵糡תΪ��˥����ѹ
    SUBROUTINE HZA(HXYZ,HHH,XYZ,WE,NS,SX)
    IMPLICIT NONE
    INTEGER WE,NS,SX,I,K1,K2,N,K3
    REAL(8) XYZ(NS*WE*SX,3),PI,U0,DRT1,TEMP,YUAN(4,3)
    REAL(8) HXYZ(WE*NS*SX,3),HHH(WE*NS*SX,3),TEMP1
    REAL(8) X(SX),Y(SX),XX(1),YY(1)
    REAL(8) X1(WE),Y1(WE),X2(NS),Y2(NS),XX1(1),XX2(1),YY1(1),YY2(1)
    HHH=0.0
    
    !HX
    DRT1=0.00001
    DO I=1,WE*NS*SX
     K2=1
    DO K1=0,SX
    IF((I-K1*WE*NS)>0)THEN
    IF(XYZ(I-K1*WE*NS,1)==XYZ(I,1).AND.XYZ((I-K1*WE*NS),2)==XYZ(I,2))THEN
    X(K2)=XYZ(I-K1*WE*NS,3)
    Y(K2)=HXYZ(I-K1*WE*NS,2)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+(K1+1)*WE*NS)<=WE*NS*SX)THEN
    IF(XYZ(I+(K1+1)*WE*NS,2)==XYZ(I,2).AND.XYZ((I+(K1+1)*WE*NS),1)==XYZ(I,1))THEN
    X(K2)=XYZ(I+(K1+1)*WE*NS,3)
    Y(K2)=HXYZ(I+(K1+1)*WE*NS,2)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO
      DO K1=1,SX
	       DO K2=K1+1,SX
	         IF (X(K1)>X(K2))THEN
	            TEMP=X(K1)
	            X(K1)=X(K2)
	            X(K2)=TEMP
	            TEMP1=Y(K1)
	            Y(K1)=Y(K2)
	            Y(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
	      XX(1)=XYZ(I,3)
    IF(XYZ(I,3)==X(1))XX(1)=XYZ(I,3)+DRT1
    IF(XYZ(I,3)==X(SX))XX(1)=XYZ(I,3)-DRT1
    CALL SPLN(X,Y,SX,XX,YY,1,SX-1,SX-2)
    
    K2=1
    DO K1=0,NS
    IF((I-K1*WE)>0)THEN
    IF(XYZ(I-K1*WE,1)==XYZ(I,1).AND.XYZ(I-K1*WE,3)==XYZ(I,3))THEN
    X2(K2)=XYZ(I-K1*WE,2)
    Y2(K2)=HXYZ(I-K1*WE,3)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+(K1+1)*WE)<=WE*NS*SX)THEN
    IF(XYZ(I+(K1+1)*WE,1)==XYZ(I,1).AND.XYZ(I+(K1+1)*WE,3)==XYZ(I,3))THEN
    X2(K2)=XYZ(I+(K1+1)*WE,2)
    Y2(K2)=HXYZ(I+(K1+1)*WE,3)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO
      DO K1=1,NS
	       DO K2=K1+1,NS
	         IF (X2(K1)>X2(K2))THEN
	            TEMP=X2(K1)
	            X2(K1)=X2(K2)
	            X2(K2)=TEMP
	            TEMP1=Y2(K1)
	            Y2(K1)=Y2(K2)
	            Y2(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
      XX2(1)=XYZ(I,2)
    IF(XYZ(I,2)==X2(1))XX2(1)=XYZ(I,2)+DRT1
    IF(XYZ(I,2)==X2(NS))XX2(1)=XYZ(I,2)-DRT1
    CALL SPLN(X2,Y2,NS,XX2,YY2,1,NS-1,NS-2)
    
    HHH(I,1)=YY(1)-YY2(1)
    ENDDO
    !HY
    DO I=1,WE*NS*SX
     K2=1
    DO K1=0,SX
    IF((I-K1*WE*NS)>0)THEN
    IF(XYZ(I-K1*WE*NS,1)==XYZ(I,1).AND.XYZ((I-K1*WE*NS),2)==XYZ(I,2))THEN
    X(K2)=XYZ(I-K1*WE*NS,3)
    Y(K2)=HXYZ(I-K1*WE*NS,1)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+(K1+1)*WE*NS)<=WE*NS*SX)THEN
    IF(XYZ(I+(K1+1)*WE*NS,2)==XYZ(I,2).AND.XYZ((I+(K1+1)*WE*NS),1)==XYZ(I,1))THEN
    X(K2)=XYZ(I+(K1+1)*WE*NS,3)
    Y(K2)=HXYZ(I+(K1+1)*WE*NS,1)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO
      DO K1=1,SX
	       DO K2=K1+1,SX
	         IF (X(K1)>X(K2))THEN
	            TEMP=X(K1)
	            X(K1)=X(K2)
	            X(K2)=TEMP
	            TEMP1=Y(K1)
	            Y(K1)=Y(K2)
	            Y(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
	      XX(1)=XYZ(I,3)
    IF(XYZ(I,3)==X(1))XX(1)=XYZ(I,3)+DRT1
    IF(XYZ(I,3)==X(SX))XX(1)=XYZ(I,3)-DRT1
    CALL SPLN(X,Y,SX,XX,YY,1,SX-1,SX-2) 
     K2=1
    DO K1=0,WE
    IF((I-K1)>0)THEN
    IF(XYZ(I-K1,2)==XYZ(I,2).AND.XYZ((I-K1),3)==XYZ(I,3))THEN
    X1(K2)=XYZ(I-K1,1)
    Y1(K2)=HXYZ(I-K1,3)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+K1+1)<=WE*NS*SX)THEN
    IF(XYZ(I+K1+1,2)==XYZ(I,2).AND.XYZ((I+K1+1),3)==XYZ(I,3))THEN
    X1(K2)=XYZ(I+K1+1,1)
    Y1(K2)=HXYZ(I+K1+1,3)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO 
      DO K1=1,WE
	       DO K2=K1+1,WE
	         IF (X1(K1)>X1(K2))THEN
	            TEMP=X1(K1)
	            X1(K1)=X1(K2)
	            X1(K2)=TEMP
	            TEMP1=Y1(K1)
	            Y1(K1)=Y1(K2)
	            Y1(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
      XX1(1)=XYZ(I,1)
    IF(XYZ(I,1)==X1(1))XX1(1)=XYZ(I,1)+DRT1
    IF(XYZ(I,1)==X1(WE))XX1(1)=XYZ(I,1)-DRT1
    CALL SPLN(X1,Y1,WE,XX1,YY1,1,WE-1,WE-2)
    HHH(I,2)=YY1(1)-YY(1)  
    ENDDO
    !HZ
    DO I=1,WE*NS*SX 
    K2=1
    DO K1=0,WE
    IF((I-K1)>0)THEN
    IF(XYZ(I-K1,2)==XYZ(I,2).AND.XYZ((I-K1),3)==XYZ(I,3))THEN
    X1(K2)=XYZ(I-K1,1)
    Y1(K2)=HXYZ(I-K1,2)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+K1+1)<=WE*NS*SX)THEN
    IF(XYZ(I+K1+1,2)==XYZ(I,2).AND.XYZ((I+K1+1),3)==XYZ(I,3))THEN
    X1(K2)=XYZ(I+K1+1,1)
    Y1(K2)=HXYZ(I+K1+1,2)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO 
      DO K1=1,WE
	       DO K2=K1+1,WE
	         IF (X1(K1)>X1(K2))THEN
	            TEMP=X1(K1)
	            X1(K1)=X1(K2)
	            X1(K2)=TEMP
	            TEMP1=Y1(K1)
	            Y1(K1)=Y1(K2)
	            Y1(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
      XX1(1)=XYZ(I,1)
    IF(XYZ(I,1)==X1(1))XX1(1)=XYZ(I,1)+DRT1
    IF(XYZ(I,1)==X1(WE))XX1(1)=XYZ(I,1)-DRT1
    CALL SPLN(X1,Y1,WE,XX1,YY1,1,WE-1,WE-2)
    K2=1
    DO K1=0,NS
    IF((I-K1*WE)>0)THEN
    IF(XYZ(I-K1*WE,1)==XYZ(I,1).AND.XYZ(I-K1*WE,3)==XYZ(I,3))THEN
    X2(K2)=XYZ(I-K1*WE,2)
    Y2(K2)=HXYZ(I-K1*WE,1)
    K2=K2+1
    ENDIF
    ENDIF
    IF((I+(K1+1)*WE)<=WE*NS*SX)THEN
    IF(XYZ(I+(K1+1)*WE,1)==XYZ(I,1).AND.XYZ(I+(K1+1)*WE,3)==XYZ(I,3))THEN
    X2(K2)=XYZ(I+(K1+1)*WE,2)
    Y2(K2)=HXYZ(I+(K1+1)*WE,1)
    K2=K2+1
    ENDIF
    ENDIF
    ENDDO
      DO K1=1,NS
	       DO K2=K1+1,NS
	         IF (X2(K1)>X2(K2))THEN
	            TEMP=X2(K1)
	            X2(K1)=X2(K2)
	            X2(K2)=TEMP
	            TEMP1=Y2(K1)
	            Y2(K1)=Y2(K2)
	            Y2(K2)=TEMP1
	          ENDIF
	       ENDDO
	     ENDDO
      XX2(1)=XYZ(I,2)
    IF(XYZ(I,2)==X2(1))XX2(1)=XYZ(I,2)+DRT1
    IF(XYZ(I,2)==X2(NS))XX2(1)=XYZ(I,2)-DRT1
    CALL SPLN(X2,Y2,NS,XX2,YY2,1,NS-1,NS-2)
    HHH(I,3)=YY2(1)-YY1(1)
    ENDDO
    
    ENDSUBROUTINE
!*******************************************************************
!                 ����������ֵ����
!         n-----��ֵ�ڵ����(n>2);
!         m-----����ֵ�������ĵ�����
!         x(n)--��Ÿ����Ĳ�ֵ�ڵ㣬Ҫ��x(i)<x(i+1),i=1,2,...n-1��
!         y(n)--��Ÿ����ڵ㴦�ĺ���ֵ��
!         t(m)--��Ų�ֵ�㣬Ҫ��x(1)<t(1)<t(2)<....<t(m)<x(n).
!         f(m)--��Ų�ֵ�����
!         f1(m)--��Ų�ֵ��һ�׵���f'(m).
!*******************************************************************
    SUBROUTINE SPLN(X,Y,N,T,F1,M,N1,N2)
    IMPLICIT NONE
    INTEGER N,M,N1,N2,I,K,J
    REAL(8) X(N),T(M)
	REAL(8) Y(N),F(M),F1(M),S2(N),H(N1),DY(N1),S(N1)
	REAL(8) E(N2),Z,H1,H2,H3,H4
	DO  I=1,N1
	H(I)=X(I+1)-X(I)
	DY(I)=(Y(I+1)-Y(I))/H(I)
    ENDDO
	S2(1)=0.
	S2(N)=0.
	DO I=2,N1
 	S2(I)=6.*(DY(I)-DY(I-1))
    ENDDO
	Z=0.5/(H(1)+H(2))
	S(1)=-H(2)*Z
	E(1)=S2(2)*Z
	DO I=2,N2
	K=I-1
	J=I+1
	Z=1./(2.*(H(I)+H(J))+H(I)*S(K))
	S(I)=-H(J)*Z
	E(I)=(S2(J)-H(I)*E(K))*Z
    ENDDO
	S2(N1)=E(N2)
	DO I=N2,2,-1
	K=I-1
	S2(I)=S(K)*S2(I+1)+E(K)
    ENDDO
	DO I=1,N1
 	S(I)=(S2(I+1)-S2(I))/H(I)
    ENDDO
	I=2
	K=1
	DO J=1,M
100	IF(T(J).GT.X(I)) then
	K=I
	I=I+1
	GOTO 100
	ELSE
	H1=T(J)-X(K)
	H2=T(J)-X(I)
	H3=H1*H2
	H4=S2(k)+H1*S(K)
	Z=(S2(I)+S2(K)+H4)/6.
	F(J)=Y(K)+H1*DY(K)+H3*Z
	F1(J)=DY(K)+Z*(H1+H2)+H3*S(K)/6.
	ENDIF
    ENDDO
	RETURN
	ENDSUBROUTINE

