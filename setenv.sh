# Get the path of scrip dir
SCRIPT_DIR=$(cd $(dirname $0); pwd)

export PATH=$PATH:$(SCRIPT_DIR)/bins
export F3DS_LIBS=$(SCRIPT_DIR)/libs
export F3DS_MODS=$(SCRIPT_DIR)/mods