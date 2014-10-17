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

DEBUG = -02 -Wall

SELF = basic.mk

CFLAGS = $(DEBUG)

INC_FLAGS = -I./
# INC_FLAGS += -I$(HOME)/include

.c.o :
	$(CC) $(CFLAGS) $(INC_FLAGS) -c $<

SDE_WA_SRC = sde_wa.c
SDE_WA_SRC += sde_wa_em.c
SDE_WA_SRC += sde_wa_nv.c
SDE_WA_SRC += sde_wa_nn.c
SDE_WA_SRC += sde_wa_butcher.c
SDE_WA_OBJS = ${SDE_WA_SRC:.c=.o}

SDE_WA = sde_wa
TARGET_SDE_WA = lib$(SDE_WA).a
TARGETS = $(TARGET_SDE_WA) 

$(TARGET_SDE_WA): $(SDE_WA_OBJS)
	rm -f $@
	ar cvr $@ $(SDE_WA_OBJS)
clean :
	rm -f *.o $(TARGETS)

all: $(TARGETS)

Makefile : $(SELF)
	rm -f $@
	cp $(SELF) $@
	chmod +w $@
	echo '# Automatically-generated dependencies list:' >>$@
	gcc ${CFLAGS} -MM	\
	${SDE_WA_SRC}	\
	>> $@
	chmod -w $@
# Automatically-generated dependencies list:
sde_wa.o: sde_wa.c sde_wa.h sde_wa_errno.h sde_wa_em.h sde_wa_nv.h \
  sde_wa_nn.h
sde_wa_em.o: sde_wa_em.c sde_wa.h sde_wa_errno.h sde_wa_em.h
sde_wa_nv.o: sde_wa_nv.c sde_wa.h sde_wa_errno.h sde_wa_nv.h \
  sde_wa_butcher.h
sde_wa_nn.o: sde_wa_nn.c sde_wa.h sde_wa_errno.h sde_wa_nn.h \
  sde_wa_butcher.h
sde_wa_butcher.o: sde_wa_butcher.c sde_wa_butcher.h sde_wa.h \
  sde_wa_errno.h sde_wa_nn.h
