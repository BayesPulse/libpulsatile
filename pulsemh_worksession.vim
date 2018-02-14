let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +79 include/singlesubject/ss_drawlocations.h
badd +0 tests/ss_draws_test.cpp
badd +1 tests/proposalvariance_tests.cpp
badd +0 include/singlesubject/ss_drawfixedeffects.h
badd +0 tests/mh_tests.cpp
badd +4 include/mh.h
badd +0 term://.//31347:/bin/bash
badd +0 include/utils.h
badd +0 include/datastructures.h
badd +0 include/Pati
badd +0 include/patient.h
argglobal
silent! argdel *
$argadd include/singlesubject/ss_drawlocations.h
set stal=2
edit tests/ss_draws_test.cpp
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
exe 'vert 1resize ' . ((&columns * 31 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 165 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 165 + 181) / 363)
argglobal
enew
file NERD_tree_1
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
let s:l = 14 - ((13 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
14
normal! 035|
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
let s:l = 85 - ((39 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
85
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 165 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 165 + 181) / 363)
tabedit include/singlesubject/ss_drawlocations.h
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
let s:l = 80 - ((78 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
80
normal! 0
wincmd w
argglobal
if bufexists('tests/mh_tests.cpp') | buffer tests/mh_tests.cpp | else | edit tests/mh_tests.cpp | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 261 - ((83 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
261
normal! 0
wincmd w
argglobal
if bufexists('include/singlesubject/ss_drawfixedeffects.h') | buffer include/singlesubject/ss_drawfixedeffects.h | else | edit include/singlesubject/ss_drawfixedeffects.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 25 - ((24 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
25
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 121 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 120 + 181) / 363)
exe 'vert 3resize ' . ((&columns * 120 + 181) / 363)
tabedit include/mh.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 181 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 181 + 181) / 363)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 108 - ((74 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
108
normal! 038|
wincmd w
argglobal
if bufexists('term://.//31347:/bin/bash') | buffer term://.//31347:/bin/bash | else | edit term://.//31347:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 1016 - ((83 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1016
normal! 037|
wincmd w
exe 'vert 1resize ' . ((&columns * 181 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 181 + 181) / 363)
tabedit include/datastructures.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 181 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 181 + 181) / 363)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 22 - ((21 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
22
normal! 02|
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
let s:l = 20 - ((19 * winheight(0) + 42) / 84)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
20
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 181 + 181) / 363)
exe 'vert 2resize ' . ((&columns * 181 + 181) / 363)
tabnext 2
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
