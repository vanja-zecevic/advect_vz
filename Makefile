#  Copyright (C) 2012 Vanja Zecevic
#  Contact vanja.zecevic@sydney.uni.edu.au
#
#  This file is part of advect_vz
#
#  advect_vz is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  advect_vz is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 

NAME=advection

CC = gcc
CFLAGS = -g -Wall -O3 -fopenmp -lm -lrt -Winline --param\
 max-inline-insns-single=5000 --param inline-unit-growth=500
INCLUDES = -I./ -I./src/

IN_ADV := \
  flux.c \
  core.c \
  core_test.c \
  tools.c \
  bc.c \
  init.c

IN_LIBVZ := \
  exit_vz.c \
  lliffe_vz.c \
  cfg_vz.c

IN_RUN := \
  run2d.c

GEN_SRC := \
  src/adv/core.c \
  src/adv/flux.c 
GEN_HEADERS := $(GEN_SRC:%.c=%.h)

OBJ_ADV    := $(IN_ADV:%.c=obj/adv/%.o)
OBJ_LIBVZ  := $(IN_LIBVZ:%.c=obj/libvz/%.o)
OBJ_RUN    := $(IN_RUN:%.c=obj/run/%.o)

DEP_ADV    := $(IN_ADV:%.c=dep/adv/%.d)
DEP_LIBVZ  := $(IN_LIBVZ:%.c=dep/libvz/%.d)
DEP_RUN    := $(IN_RUN:%.c=dep/run/%.d)

#-----------------------------------------------------------------------------
all: $(IN_RUN:%.c=bin/%) bin/gen_vz

bin/%: $(OBJ_ADV) $(OBJ_LIBVZ) $(OBJ_RUN)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@

#-----------------------------------------------------------------------------
# Objects.
obj/adv/%.o: src/adv/%.c
	@echo -e '\033[1;34mCC CFLAGS INCLUDES\033[00m' -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

obj/libvz/%.o: src/libvz/%.c
	@echo -e '\033[1;34mCC CFLAGS INCLUDES\033[00m' -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

obj/run/%.o: src/run/%.c
	@echo -e '\033[1;34mCC CFLAGS INCLUDES\033[00m' -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

#-----------------------------------------------------------------------------
# Dependencies.
dep/adv/%.d: src/adv/%.c $(GEN_HEADERS)
	@echo -e '\033[1;36mCC CFLAGS INCLUDES -MM -MT\033[00m' $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -MM -MT $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<

dep/libvz/%.d: src/libvz/%.c $(GEN_HEADERS)
	@echo -e '\033[1;36mCC CFLAGS INCLUDES -MM -MT\033[00m' $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -MM -MT $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<

dep/run/%.d: src/run/%.c  $(GEN_HEADERS)
	@echo -e '\033[1;36mCC CFLAGS INCLUDES -MM -MT\033[00m' $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<
	@$(CC) $(CFLAGS) $(INCLUDES) -MM -MT $@ \
  -MT $(subst dep/, obj/, $(@:.d=.o)) -o $@ $<

#-----------------------------------------------------------------------------
# Generated files.
src/adv/%.c: src/adv/%.gen.c bin/gen_vz
	cat $< | bin/gen_vz > $@

src/adv/%.h: src/adv/%.gen.h bin/gen_vz
	cat $< | bin/gen_vz > $@

#-----------------------------------------------------------------------------
# gen_vz
bin/gen_vz: src/libvz/gen_vz.c
	gcc -Wall $< -o $@ -lfl

src/libvz/gen_vz.c: src/libvz/gen_vz.lex
	flex -s -o$@ $<
#-----------------------------------------------------------------------------
# Cleaning.
clean:
	rm -vf bin/*
	find obj -type f -print0 | xargs -0 rm -vf
	find dep -type f -print0 | xargs -0 rm -vf
	rm -vf src/libvz/gen_vz.c
	rm -vf $(GEN_SRC) $(GEN_HEADERS)

#-----------------------------------------------------------------------------

ifneq ($(MAKECMDGOALS),clean)
  include $(DEP_ADV)
  include $(DEP_LIBVZ)
  include $(DEP_RUN)
endif

