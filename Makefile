HOME = $(shell echo $$HOME)

SRC = /fastdata/bo1mesx/workshop/
BIN = $(HOME)genabel/

.PHONY: all
all:  $(DIR)phe_RUFF.txt $(DIR)gen_RUFF_qc.raw

$(BIN)%: $(SRC)%
  mkdir -p $(BIN)
  cp $^ $@
