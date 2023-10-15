      REAL(10) function dcabs1(z)
      double complex z,zz
      REAL(10) t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = abs(t(1)) + abs(t(2))
      return
      end
