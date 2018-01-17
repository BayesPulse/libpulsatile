let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +36 include/mh.h
badd +2 include/counter.h
badd +56 include/datastructures.h
badd +1 include/patient.h
badd +1 include/population.h
badd +105 include/proposalvariance.h
badd +8 include/utils.h
badd +21 src/tempmain.cpp
badd +3 tests/datastructures_tests.cpp
badd +0 ../Software/pulsatile/src/mcmc.c
badd +2 src/ss_random_effects.cpp
badd +72 tmp-draw_random_effects.c
badd +1 tmp-draw_fixed_effects.cpp
badd +54 tmp-draw_fixed_effects.c
badd +72 src/ss_fixed_effects.cpp
badd +14 devtmp/parm_location_map.R
badd +0 tests/ss_draws_test.cpp
badd +0 tests/patient_tests.cpp
badd +44 include/ss_fixed_effects.h
badd +0 tests/proposalvariance_tests.cpp
argglobal
silent! argdel *
$argadd include/mh.h
set stal=2
edit include/mh.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 214) / 428)
exe '2resize ' . ((&lines * 47 + 48) / 96)
exe 'vert 2resize ' . ((&columns * 110 + 214) / 428)
exe '3resize ' . ((&lines * 47 + 48) / 96)
exe 'vert 3resize ' . ((&columns * 142 + 214) / 428)
exe '4resize ' . ((&lines * 45 + 48) / 96)
exe 'vert 4resize ' . ((&columns * 126 + 214) / 428)
exe '5resize ' . ((&lines * 45 + 48) / 96)
exe 'vert 5resize ' . ((&columns * 126 + 214) / 428)
exe 'vert 6resize ' . ((&columns * 142 + 214) / 428)
argglobal
enew
file NERD_tree_2
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
let s:l = 46 - ((31 * winheight(0) + 23) / 47)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
46
normal! 0
wincmd w
argglobal
if bufexists('include/counter.h') | buffer include/counter.h | else | edit include/counter.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 29 - ((28 * winheight(0) + 23) / 47)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
29
normal! 038|
wincmd w
argglobal
if bufexists('include/proposalvariance.h') | buffer include/proposalvariance.h | else | edit include/proposalvariance.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=1
setlocal fen
29
normal! zo
119
normal! zo
let s:l = 59 - ((39 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
59
normal! 07|
wincmd w
argglobal
if bufexists('include/ss_fixed_effects.h') | buffer include/ss_fixed_effects.h | else | edit include/ss_fixed_effects.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 21 - ((13 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
21
normal! 0
wincmd w
argglobal
if bufexists('tests/ss_draws_test.cpp') | buffer tests/ss_draws_test.cpp | else | edit tests/ss_draws_test.cpp | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 99 - ((87 * winheight(0) + 46) / 93)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
99
normal! 025|
wincmd w
5wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 214) / 428)
exe '2resize ' . ((&lines * 47 + 48) / 96)
exe 'vert 2resize ' . ((&columns * 110 + 214) / 428)
exe '3resize ' . ((&lines * 47 + 48) / 96)
exe 'vert 3resize ' . ((&columns * 142 + 214) / 428)
exe '4resize ' . ((&lines * 45 + 48) / 96)
exe 'vert 4resize ' . ((&columns * 126 + 214) / 428)
exe '5resize ' . ((&lines * 45 + 48) / 96)
exe 'vert 5resize ' . ((&columns * 126 + 214) / 428)
exe 'vert 6resize ' . ((&columns * 142 + 214) / 428)
tabedit tests/proposalvariance_tests.cpp
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
let s:l = 66 - ((18 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
66
normal! 0
tabedit src/ss_random_effects.cpp
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
let s:l = 3 - ((2 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
3
normal! 0
tabedit src/ss_fixed_effects.cpp
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
let s:l = 2 - ((1 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
tabedit include/ss_fixed_effects.h
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
let s:l = 2 - ((1 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
tabedit tests/datastructures_tests.cpp
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
let s:l = 448 - ((67 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
448
normal! 017|
tabedit src/ss_random_effects.cpp
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
let s:l = 23 - ((22 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
23
normal! 0
tabedit include/datastructures.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 213 + 214) / 428)
exe 'vert 2resize ' . ((&columns * 214 + 214) / 428)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 235 - ((59 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
235
normal! 05|
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
let s:l = 72 - ((71 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
72
normal! 03|
wincmd w
exe 'vert 1resize ' . ((&columns * 213 + 214) / 428)
exe 'vert 2resize ' . ((&columns * 214 + 214) / 428)
tabedit ../Software/pulsatile/src/mcmc.c
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
let s:l = 1093 - ((70 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1093
normal! 0
tabedit include/proposalvariance.h
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 6 - ((5 * winheight(0) + 47) / 95)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
6
normal! 0
tabnext 1
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
