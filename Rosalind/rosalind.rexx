#!/usr/bin/rexx

/* call HEA 'hea.txt' */
/* call HS 'hs.txt' */
call PS 'ps.txt'

exit 0

PS: procedure

  parse arg f

  /* Read data */
  arr.0 = linein(f)
  call lineToInts linein(f)
  k = linein(f)

  do i = 1 to k

    /* Find smallest value */
    m = i
    do j = i + 1 to arr.0
      if arr.j < arr.m then
        m = j
    end

    /* Swap */
    tmp = arr.i
    arr.i = arr.m
    arr.m = tmp

  end

  /* Show result */
  do i = 1 to k
    call charout , arr.i' '
  end
  say ''

  return

HS: procedure expose arr.

  parse arg f

  /* Create heap */
  call HEA f

  hSize = arr.0
  do i = arr.0 to 2 by -1

    /* Swap last element with root */
    tmp = arr.hSize
    arr.hSize = arr.1
    arr.1 = tmp

    /* Decrease heap size */
    hSize = hSize - 1

    cIndex = 1
    do forever

      l = cIndex * 2
      r = l + 1

      if l <= hSize & arr.l > arr.cIndex & (r > hSize | arr.r <= arr.l) then

        /* Move left */
        do
          tmp = arr.cIndex
          arr.cIndex = arr.l
          arr.l = tmp
          cIndex = l
        end

      else if r <= hSize & arr.r > arr.cIndex then

        /* Move right */
        do
          tmp = arr.cIndex
          arr.cIndex = arr.r
          arr.r = tmp
          cIndex = r
        end

      else

        leave
      
    end

  end

  /* Result */
  /*
  do i = 1 to arr.0
    call charout , arr.i' '
  end
  say ''
  */

  return

HEA: procedure expose arr.

  parse arg f

  /* Set array size */
  arr.0 = linein(f)

  /* Read array from file */
  call lineToInts linein(f)

  /* Close file */
  call lineout f

  /* Create heap */
  do i = 2 to arr.0

    cIndex = i
    do while cIndex > 1

      cVal = arr.cIndex
      parIndex = cIndex % 2
      parVal = arr.parIndex

      if parVal < cVal then do
        arr.cIndex = parVal
        arr.parIndex = cVal
      end

      cIndex = parIndex

    end

  end

  /* Display result */
  /*
  do i = 1 to arr.0
    call charout , arr.i' '
  end
  say ''
  */

  return

lineToInts: procedure expose arr.

  arg line

  cString = '' ; arrIndex = 1

  ll = length(line)
  do i = 1 to ll

    s = substr(line,i,1)

    if ' ' \= s then
      cString = cString || s
    else do
      arr.arrIndex = cString
      arrIndex = arrIndex + 1
      cString = ''
    end

  end

  arr.arrIndex = cString

  return
