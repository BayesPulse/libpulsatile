let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +70 include/singlesubject/ss_drawlocations.h
badd +1 tests/ss_draws_test.cpp
badd +1 tests/proposalvariance_tests.cpp
badd +41 include/singlesubject/ss_drawfixedeffects.h
badd +1 tests/mh_tests.cpp
badd +85 include/mh.h
badd +97 term://.//28642:/bin/bash
badd +86 include/utils.h
badd +1 include/datastructures.h
badd +1 include/Pati
badd +1 include/patient.h
badd +3 include/population.h
badd +26 ~/Projects/BayesPulse/pulsatile/src/mcmc.c
badd +1 ~/Projects/BayesPulse/pulsatile/R/mcmc_pulseparameters.R
badd +83 ~/Projects/BayesPulse/pulsatile/R/pulsespec.R
badd +73 include/proposalvariance.h
badd +36 term://.//26262:/bin/bash
badd +3 test.out
badd +0 term://.//19608:/bin/bash
argglobal
silent! argdel *
$argadd include/singlesubject/ss_drawlocations.h
set stal=2
edit tests/mh_tests.cpp
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
exe 'vert 1resize ' . ((&columns * 121 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 120 + 181) / 363)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 280 - ((80 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
280
normal! 0
wincmd w
argglobal
if bufexists('include/singlesubject/ss_drawlocations.h') | buffer include/singlesubject/ss_drawlocations.h | else | edit include/singlesubject/ss_drawlocations.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 42 - ((41 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
42
normal! 0
wincmd w
argglobal
if bufexists('include/mh.h') | buffer include/mh.h | else | edit include/mh.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 80 - ((79 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
80
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 121 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 120 + 181) / 363)
tabedit test.out
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
exe 'vert 1resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 121 + 181) / 363)
argglobal
if bufexists('term://.//26262:/bin/bash') | buffer term://.//26262:/bin/bash | else | edit term://.//26262:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 51 - ((50 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
51
normal! 0
wincmd w
argglobal
if bufexists('term://.//19608:/bin/bash') | buffer term://.//19608:/bin/bash | else | edit term://.//19608:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 10084 - ((83 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10084
normal! 037|
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
let s:l = 2 - ((1 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 121 + 181) / 363)
tabedit include/datastructures.h
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
exe 'vert 1resize ' . ((&columns * 175 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 1 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 185 + 181) / 363)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 224 - ((54 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
224
normal! 03|
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
let s:l = 188 - ((0 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
188
normal! 0
wincmd w
argglobal
if bufexists('include/population.h') | buffer include/population.h | else | edit include/population.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 72 - ((71 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
72
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 175 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 1 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 185 + 181) / 363)
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
