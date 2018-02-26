let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +76 include/singlesubject/ss_draw_baselinehalflife.h
badd +114 include/mh.h
badd +23 include/proposalvariance.h
badd +6 include/singlesubject/ss_draw_fixedeffects.h
badd +162 src/tempmain.cpp
badd +448 tests/mh_tests.cpp
badd +75 include/singlesubject/ss_draw_locations.h
badd +142 tests/proposalvariance_tests.cpp
badd +68 include/utils.h
badd +100 include/datastructures.h
badd +1391 ~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c
badd +3 include/patient.h
badd +30 include/singlesubject/ss_draw_error.h
badd +92 include/singlesubject/ss_draw_randomeffects.h
badd +0 term://.//52974:/bin/bash
badd +269 ~/Projects/BayesPulse/Software/pulsatile/src/birth_death.c
badd +0 include/birthdeath.h
badd +35 ~/Projects/BayesPulse/Software/jointpulsatile/JointModelCode_restructured/src/bd_trigger.c
badd +153 tests/patient_tests.cpp
badd +0 include/singlesubject/ss_draw_tvarscale.h
badd +0 ~/Projects/BayesPulse/Software/pulsatile/R/pulsespec.R
badd +0 ~/Projects/BayesPulse/Software/pulsatile/src/calculations.c
argglobal
silent! argdel *
$argadd include/singlesubject/ss_draw_baselinehalflife.h
set stal=2
edit include/singlesubject/ss_draw_locations.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 102 + 119) / 238)
exe 'vert 3resize ' . ((&columns * 103 + 119) / 238)
argglobal
enew
file NERD_tree_6
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal nofen
wincmd w
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 99 - ((50 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
99
normal! 0
wincmd w
argglobal
if bufexists('include/patient.h') | buffer include/patient.h | else | edit include/patient.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 43 - ((12 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
43
normal! 03|
wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 102 + 119) / 238)
exe 'vert 3resize ' . ((&columns * 103 + 119) / 238)
tabedit include/birthdeath.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 119 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 118 + 119) / 238)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 211 - ((50 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
211
normal! 07|
wincmd w
argglobal
if bufexists('include/utils.h') | buffer include/utils.h | else | edit include/utils.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 68 - ((33 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
68
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 119 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 118 + 119) / 238)
tabedit ~/Projects/BayesPulse/Software/pulsatile/src/calculations.c
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
if bufexists('~/Projects/BayesPulse/Software/pulsatile/src/birth_death.c') | buffer ~/Projects/BayesPulse/Software/pulsatile/src/birth_death.c | else | edit ~/Projects/BayesPulse/Software/pulsatile/src/birth_death.c | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=1
setlocal fen
82
normal! zo
498
normal! zo
648
normal! zo
712
normal! zo
let s:l = 286 - ((45 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
286
normal! 016|
wincmd w
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
tabedit include/singlesubject/ss_draw_tvarscale.h
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 35 - ((34 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
35
normal! 0
tabedit include/patient.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=1
setlocal fen
32
normal! zo
let s:l = 92 - ((35 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
92
normal! 0
wincmd w
argglobal
if bufexists('include/datastructures.h') | buffer include/datastructures.h | else | edit include/datastructures.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 183 - ((17 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
183
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
tabedit ~/Projects/BayesPulse/Software/pulsatile/R/pulsespec.R
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 157 - ((27 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
157
normal! 062|
tabedit include/mh.h
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 46 - ((17 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
46
normal! 03|
tabedit include/patient.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
argglobal
if bufexists('term://.//52974:/bin/bash') | buffer term://.//52974:/bin/bash | else | edit term://.//52974:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 713 - ((55 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
713
normal! 038|
wincmd w
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 150 - ((29 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
150
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 118 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 238)
tabnext 8
set stal=1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
