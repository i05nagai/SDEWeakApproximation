#
# makefile for libsde_wa.a
#
# $Source$
# $Author$
# $Date$
# $Revision$
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; version 2.1 of the License.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License version 2.1 for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# version 2.1 along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, US
#
.SUFFIXES :
.SUFFIXES : .o .c


CC = gcc
DEBUG = -O2 -Wall
CFLAGS = $(DEBUG)

SELF = basic.mk
MAKEFILE = Makefile
DIRS = $(shell find . -mindepth 1 -maxdepth 1 -type d | grep -v "\/\.")
SELFS = $(addsuffix /$(SELF), $(DIRS))
MAKEFILES = $(addsuffix /$(MAKEFILE), $(DIRS))

SDE_WA_SELFS = $(shell find . -name '$(SELF)')
SDE_WA_MAKEFILES = $(shell find . -name '$(MAKEFILE)')
SDE_WA_SRC = $(shell find . -name '*.c')
SDE_WA_OBJS =  $(SDE_WA_SRC:.c=.o)
SDE_WA = sde_wa
TARGET_SDE_WA = lib$(SDE_WA).a
TARGETS = $(TARGET_SDE_WA) 


$(TARGET_SDE_WA): create-objs
	rm -f $@
	ar cvr $@ $(SDE_WA_OBJS)

Makefiles: $(SELFS)
	for dir in $(DIRS) ; do	\
	 cd $$dir	;\
	 make -f $(SELF) $(MAKEFILE)	;\
	 cd ..	;\
	done

create-objs: $(SDE_WA_SELFS)
	for dir in $(DIRS); do	\
		cd $$dir	;\
		make	;\
		cd ..	;\
	done


.PHONY: all
all: clean $(TARGETS)

.PHONY: clean
clean :
	rm -f $(SDE_WA_OBJS) $(TARGETS)

