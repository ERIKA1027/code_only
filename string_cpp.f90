module string
  implicit none
  private
  integer,parameter :: CHARLEN = 1000 
  interface num2char
     module procedure f2char, i2char, ill2char
  end interface num2char
   public :: CHARLEN, num2char, i2char0, concat, strcmprs
  contains
    function f2char(num, charformat) result(char)
      character(len=CHARLEN) :: char
      real(kind=8),intent(in) :: num
      character(len=*),intent(IN),optional :: charformat
      character(len=100) :: chformat
      call clean_char(char)
      if (present(charformat)) then
         chformat = charformat
      else
         chformat = '(1PE14.7)'
      endif
      write(char,chformat) num
      char = adjustl(char)
    end function f2char
    function i2char(num, charformat) result(char)
      character(len=CHARLEN) :: char
      integer,intent(in) :: num
      character(len=*),intent(IN),optional :: charformat
      character(len=100) :: chformat
      call clean_char(char)
      if (present(charformat)) then
         chformat = charformat
      else
         chformat = '(I0)'
      endif
      write(char,chformat) num
      char = adjustl(char)
    end function i2char
    function ill2char(num, charformat) result(char)
      character(len=CHARLEN) :: char
      integer(kind=8),intent(in) :: num
      character(len=*),intent(IN),optional :: charformat
      character(len=100) :: chformat
      call clean_char(char)
      if (present(charformat)) then
         chformat = charformat
      else
         chformat = '(I0)'
      endif
      write(char,chformat) num
      char = adjustl(char)
    end function ill2char
    function i2char0(num,charformat) result(char)
      integer,intent(in) :: num
      character(len=CHARLEN) :: char
      character(len=*),intent(IN),optional :: charformat
      character(len=100) :: chformat
      integer :: i
      call clean_char(char)
      if (present(charformat)) then
         chformat = charformat
      else
         chformat = '(I0)'
      endif
      write(char,chformat) num
      char = adjustr(char)
      do i=1,len(char)
         if(char(i:i).eq.' ') then
            char(i:i)='0'
         else
            return
         endif
      enddo
    end function i2char0
    subroutine clean_char(char)
      character(len=*),intent(INOUT) :: char
      integer :: n
      do n=1,len(char)
         char(n:n)=' '
      enddo
    end subroutine clean_char
    function concat(ch1,ch2) result(ch)
      character(len=*),intent(IN) :: ch1, ch2
      character(len=CHARLEN) :: ch
      integer :: len
      call clean_char(ch)
      ch = trim(ch1) // trim(ch2)
      len = len_trim(ch)
    end function concat
    function strcmprs(chin) result(chout)
      character(len=*),intent(in) :: chin
      character(len=CHARLEN) :: chout
      integer :: lenc,ipoint,i
      call clean_char(chout)
      lenc = len(chout)
      ipoint = 1
      do i=1,lenc
         if (chin(i:i) .ne. ' ') then
            chout(ipoint:ipoint)=chin(i:i)
            ipoint = ipoint + 1
         endif
      enddo
    end function strcmprs
end module string
