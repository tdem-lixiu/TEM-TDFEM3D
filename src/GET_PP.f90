  
    
    !���ӳ���õ�ÿ����Ԫ�ĵ�����
   SUBROUTINE  GET_PP(PP,NS,SX,WE)
   IMPLICIT NONE
   INTEGER NS,SX,WE,LENGTHNUM,I,J,K,N
   PARAMETER(LENGTHNUM=80)
   CHARACTER(LEN=lengthNum)::NUM       
   INTEGER LOCA1NUM,LOCA2NUM,LOCA1,LOCA2
   REAL(8) PP((NS-1)*(WE-1)*(SX-1))
   CHARACTER(80) FILENAME
   DO I=1,22*(WE-1)*(NS-1)
   PP(I)=0.00001                  !����������
   ENDDO
   DO I=23,SX-1
   WRITE(NUM,'(I2.2)') I
   LOCA1NUM=LOCA1(NUM,LENGTHNUM)   
   !�ǿ��ַ����׸�λ������ƶ�һλ����1001��Ϊ001
   LOCA2NUM=LOCA2(NUM,LENGTHNUM,LOCA1NUM)
   FILENAME=NUM(LOCA1NUM:LOCA2NUM)//'.txt'
   OPEN(101,FILE=FILENAME)
   READ(101,*)
   DO J=1,NS-1
   READ(101,*)(PP((I-1)*(WE-1)*(NS-1)+(J-1)*(WE-1)+K),K=1,WE-1)   !��ȡÿ����Ԫ�ĵ�����
   ENDDO
   CLOSE(101)
   ENDDO
   
   ENDSUBROUTINE
   !******************************
      INTEGER FUNCTION LOCA1(STR,LENGTH)     
      CHARACTER*(*) STR
      INTEGER LOCA2,LENGTH
      INTEGER I
      DO I=1,LENGTH
      IF (STR(I:I).NE.'') THEN
      GOTO 222
      ENDIF
      ENDDO
222   LOCA1=I
      RETURN 
      ENDFUNCTION
!*************************************************
      INTEGER FUNCTION LOCA2(STR,LENGTH,LOCA1)
      CHARACTER*(*) STR
      INTEGER LOCA1 ,LENGTH
      DO I=LOCA1+1,LENGTH
      IF(STR(I:I).EQ.'') THEN
      GOTO 223
      ENDIF
      ENDDO
223   LOCA2=I-1
      RETURN
      ENDFUNCTION


