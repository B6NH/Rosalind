#!/usr/bin/rexx

call HEA

exit 0

HEA: procedure

  f = 'hea.txt'

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
  do i = 1 to arr.0
    call charout , arr.i' '
  end
  say ''

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
