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

SRCDIR=src
INCDIR=include
VPATH = $(SRCDIR)
VPATH += $(INCDIR)

CC = gcc

DEBUG = -02 -Wall

SELF = basic.mk

CFLAGS = $(DEBUG)

INC_FLAGS = -I./
INC_FLAGS = -I./$(INCDIR)
# INC_FLAGS += -I$(HOME)/include

.c.o :
	$(CC) $(CFLAGS) $(INC_FLAGS) -c $<


SDE_WA_OBJS =  sde_wa.o
SDE_WA_OBJS += sde_wa_butcher.o
SDE_WA_OBJS += sde_wa_em.o
SDE_WA_OBJS += sde_wa_nn.o
SDE_WA_OBJS += sde_wa_nv.o
SDE_WA_OBJS += sde_wa_c3.o

SDE_WA_SRC =  $(addprefix $(SRCDIR)/, $(SDE_WA_OBJS:.o=.c))

SDE_WA = sde_wa
TARGET_SDE_WA = lib$(SDE_WA).a
TARGETS = $(TARGET_SDE_WA) 


$(TARGET_SDE_WA): $(SDE_WA_OBJS)
	rm -f $@
	ar cvr $@ $(SDE_WA_OBJS)

Makefile : $(SELF)
	rm -f $@
	cp $(SELF) $@
	chmod +w $@
	echo '# Automatically-generated dependencies list:' >>$@
	gcc ${CFLAGS} ${INC_FLAGS} -MM	\
	${SDE_WA_SRC}	\
	>> $@
	chmod -w $@

.PHONY: all
all: clean $(TARGETS)

.PHONY: clean
clean :
	rm -f *.o $(TARGETS)

# Automatically-generated dependencies list:


