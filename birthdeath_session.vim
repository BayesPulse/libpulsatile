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
badd +13 src/tempmain.cpp
badd +448 tests/mh_tests.cpp
badd +1 include/singlesubject/ss_draw_locations.h
badd +142 tests/proposalvariance_tests.cpp
badd +69 include/utils.h
badd +100 include/datastructures.h
badd +1391 ~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c
badd +221 include/patient.h
badd +30 include/singlesubject/ss_draw_error.h
badd +92 include/singlesubject/ss_draw_randomeffects.h
badd +1 term://.//1521:/bin/bash
badd +269 ~/Projects/BayesPulse/Software/pulsatile/src/birth_death.c
badd +166 include/birthdeath.h
badd +35 ~/Projects/BayesPulse/Software/jointpulsatile/JointModelCode_restructured/src/bd_trigger.c
badd +153 tests/patient_tests.cpp
badd +1 include/singlesubject/ss_draw_tvarscale.h
badd +1 ~/Projects/BayesPulse/Software/pulsatile/R/pulsespec.R
badd +1 ~/Projects/BayesPulse/Software/pulsatile/src/calculations.c
badd +0 NERD_tree_6
argglobal
silent! argdel *
$argadd include/singlesubject/ss_draw_baselinehalflife.h
set stal=2
edit include/singlesubject/ss_draw_locations.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
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
let s:l = 43 - ((17 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
43
normal! 03|
wincmd w
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
tabedit include/birthdeath.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=1
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 32 - ((25 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
32
normal! 010|
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
let s:l = 3 - ((2 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
3
normal! 0
wincmd w
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
tabedit ~/Projects/BayesPulse/Software/pulsatile/src/calculations.c
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
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
let s:l = 286 - ((63 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
286
normal! 016|
wincmd w
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 118 + 182) / 364)
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
let s:l = 35 - ((24 * winheight(0) + 28) / 56)
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
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 181 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 56 + 182) / 364)
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
let s:l = 183 - ((1 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
183
normal! 0
wincmd w
exe '1resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 1resize ' . ((&columns * 181 + 182) / 364)
exe '2resize ' . ((&lines * 56 + 41) / 82)
exe 'vert 2resize ' . ((&columns * 56 + 182) / 364)
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
tabedit include/birthdeath.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe 'vert 2resize ' . ((&columns * 244 + 182) / 364)
argglobal
if bufexists('term://.//1521:/bin/bash') | buffer term://.//1521:/bin/bash | else | edit term://.//1521:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 558 - ((78 * winheight(0) + 39) / 79)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
558
normal! 0
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
let s:l = 177 - ((73 * winheight(0) + 39) / 79)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
177
normal! 09|
wincmd w
exe 'vert 1resize ' . ((&columns * 119 + 182) / 364)
exe 'vert 2resize ' . ((&columns * 244 + 182) / 364)
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
