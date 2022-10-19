test_that("restart resid", {

  lobo  <- function() {
    ini({
      lkng <- log(0.02)
      ltau <- log(c(1, 34.1, 500))
      lec50 <- log(c(1, 50, 2000))
      kmax <- 0.01
      propSd <-  c(0, 0.3)
      addSd <- c(0, 50)
    })
    model({
      kng <- exp(lkng)
      tau <- exp(ltau)
      taulast <- tau
      ec50 <- exp(lec50)
      edrug <- kmax * cp/(ec50 + cp)
      tumor(0) <- tumor0
      d/dt(transit1) <- (edrug - transit1)/tau
      d/dt(transit2) <- (transit1 - transit2)/tau
      d/dt(transit3) <- (transit2 - transit3)/tau
      d/dt(transitlast) <- transit3/tau - transitlast/taulast
      d/dt(tumor) <- kng * tumor - transitlast * tumor
      tumor ~ prop(propSd) + add(addSd)
    })
  }

  prepfit <- qs::qdeserialize(qs::base91_decode("un]\"BAAA@QRtHACAAAAAAAuWeBAABdk1kus^^d8Ah9}?=Z:alBMc4Iv(F\":C?hVBAMZRxwFfBB7IB.y6FTL)yFQA@v;KlgSH3Vn~rL/,{CP/ez~`.3>$Wj$rcy==//#}Pu?\"V(RgKtf1J3qQg/yI7*1]/WV=iegVmPs3?a\":kEMu~/*zmX#;E4`i@It`]Ouu[N]T8G3!4A.^j0<Zd;mY7D`P6Z|KTj_p?r7u[IPuyFRoUcL\"6u|(G_on6g9c{ZLJ[_gE^&47rbL(#6W{EA%hQUW!][2;k<}\"(4SB{!~>M%xmxwMU%ET6#~xvX9rH!;S53gbLTWY]Fcri\"]7\"|Z^W{xobiiTc~DLN_;.Itj(INGKCupDYxEA^!GzfHO=aDW&(I)z}0*mZD\"^b.O!QdY2rVRD;~Z*HB(]G_Dpj]*0A`]+7VoGDF@,vk>jx}tFI>MVOnZojuABN9Bt\"O~V[n6U[kn|W74&xR7CL(Skn:CA)NP`||hQ%w/i+&c8$#KxsFdb4,qI\"Fl&lLg,?$eh&s{`QxtwPWi$GX<[*<0{to@[:NAy}a=O`wedEA*Abqhz2bL2sfII3ZRJR#5q~:FeBW%/F<]`(?Q:c(qc,DZ_d.&|J(NW~Q4kz;Us(7e+Z0YGMdvf.%XRgD]FA2D10sl^KxuPXvXSm+p}ndVY!3`o}Iq+M;i~mLmr1In0~ymm]K2x9g9Ij.UkBOTriq+93<po9tNj@%%W#FzA+/MMb]k3YmXS*RdjE{?pjnr%q6}&81.2&#ni.Au{pL>#eKT9uaMmvF7L^aDwL?sj>s|}[XMQdEx(yS|vIDfwqH:YXc(EfrzuplGn3`|X=ObNnD%;3(ST3tWr^D+vDG=cjKk!^:5ZfpXK5/dqF@dW6+*lb~@\"*H_t@3rRG9w|kG1!SRghm{sUOcDQ?.gYd?c;:IC~RF6lEfd4X;I3+4Y&8.x0P[h8]Ffo)$Y=7)|08Vid5@2btxd4I[0>E7~+CCX~|{Ve0qf^aQ1M4vsQN<B6]~R{LE*vH*8[Q|b[`QO%/.,mr9d<+O/PqesyX).kR)a[.s^OsSJLRZgrwt?NA#EQf)p,+wjGm!8aRz9^SKrW5`PTsjH_EVK.eZ)T]k0Md!mRLU%+}K!{BB>[uCto@(_b&&NA0eV*R}~Tp<Ey:7|>#Kz+AL~%_Wu@J&>C__{[#o3Lmxh$%R;264SAA)O;,:*:p.s:z+=pot[NkrFE@[}VoBF=P<AI)v_tBGUZy4^}6MUL|glMtnz]K/D|qUE``3x])hpz{k$OXX}C4B%m,AwVR)[N\"d$i}eLLJ;bFzpavPm|k^Ef;$C|S+GBP;CdOQSE<>v:;Cj#YR]M@a{68?;s&4Y4Zd[D`znL|kp1n)}0+Q`(/LcQM42JIA}ZFH+9>#dQLs>JX?#GT<W+#YOo[zpESq=N.@3/bpfI)5wKpkIRY>_M&#wui|N/HWs/WyJ/CSn*BaX!:]**?VL`zPu?xA|.+kn)M[xCjLHQ@8A6HX^T4hXT7F~JQw])fR/^Y>YBD([_QsK;[iYPcbW2Z\"BA9q~@t^rOgPg95OnA01xfA00E8k)eQtyH#MzpQ$4|et<eI&o|Uhx|t19ZFY$DAw5Z<,6p7MQM*N07@b$w):77pmoP{iZ4KvtvMyKYw|?{al}MZNtX<:RTKFB?PAdg=76:~4(HH*s62mDm2K:Qw^(9,*ICI.`?HiTsl8KIYV0C!0%0S:JM6(#2A)w(y(<9\"?n`k{0V!yrmQChNo8WWL=/*eAA_+&bKol/+*eiJi*e2dIW!hi.)+RHLR9O}Z69<<eE;aD%#B0ql|30pbBs3Z^0(q$nWfGs|LD$erJR2r@>6W(eGPRmx0B?bGrM@DA31ha^UZWMe52Egm6nn3=lcs3[gq!vBvrKiTSf],3$wjk7^`BPE#{kPh.CcbM*h:9+XPRfGx1P^)X!utw{#;**|LZtNdAwx\"IC?m[]Wil9!m3Tx}{Xx@9;`=+ODLPp)w]zc\"D.rl7/d:ismtgiesu]/CF[9v7ICZ:U?wky`@wQn])q@|Dd#IfQbwAQF`7]{f;hJCK=Vv)N6M4x4>G.u%o7I5fCay=\"^,?M!600bO\"MBTL9~4}eDkY:YkGc6G9oMzo3<9[Ye`GFvPyXiw5$Z]*PByzRza9ik5LJ;}@p;JH}+.,1om5gxAkHdn%`#kbKMJuiDOa#*M!h.2dUje501>,MIW$207SwS$M~FCMaG:`&l?j<c+ZMyFt+u:CpBWOZL)$@WywmY|LCGPdVv=>K]e\"!>HG]Tw0Xbm7SOIjt.B%cH5Ewy_YgfI6(#$I]1Q4Jte1:^W%XdqR4C#OV#._I`G$k|#>NfxcZp!Y]YnX~yZ0nx),%jfE\"u7yS,:y&oJV`;4nMg0g[a]Tidjt|y[*~_8H;*X61T/1kv9x$f/3N?~>xL$mO^0Nk6W|U[wD"))


  sim <- nlmixr2(lobo, data=prepfit, "focei", control=foceiControl(print=0))
  
})
