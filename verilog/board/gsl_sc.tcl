# ============================================================
# gsl_sc.tcl — System Console commands for the GSL board (Milestone C3)
#
# In System Console (Quartus: Tools > System Debugging Tools > System Console,
# or `system-console` from a Nios V Command Shell), in verilog/board:
#
#   source gsl_sc.tcl
#   gsl_info                          ;# detector ID, version, status
#   gsl_load rec/mb_s1.bin            ;# load a recording; the Nios V program starts it
#   gsl_wait_dump log/s1.bin          ;# wait for the result log and save it
#   gsl_run s1                        ;# both: rec/mb_s1.bin -> log/s1.bin
#
# Addresses are as seen by the JTAG-to-Avalon master: SDRAM at 0x0C000000,
# detector registers at 0x08000000 (readme.md, step 2).
# ============================================================

set GSL_SDRAM 0x0C000000
set GSL_REGS  0x08000000
set GSL_REC_MAGIC 0x47534C52
set GSL_LOG_MAGIC 0x47534C4C
set GSL_LOG_OFF   0x02000000

proc gsl_master {} {
    set paths [get_service_paths master]
    set pick ""
    foreach p $paths {
        if {[string match -nocase "*jtag*master*" $p] || [string match "*phy_*/master" $p]} { set pick $p }
    }
    if {$pick eq ""} { set pick [lindex $paths 0] }
    if {$pick eq ""} { error "no JTAG master found: is the board connected and programmed?" }
    return [claim_service master $pick gsl_lib]
}

proc gsl_rd {c addr} {
    return [expr {[lindex [master_read_32 $c $addr 1] 0]}]
}

proc gsl_info {} {
    global GSL_REGS
    set c [gsl_master]
    set id  [gsl_rd $c $GSL_REGS]
    set ver [gsl_rd $c [expr {$GSL_REGS + 0x4}]]
    set st  [gsl_rd $c [expr {$GSL_REGS + 0xC}]]
    close_service master $c
    puts [format "gsl_mm ID 0x%08X (%s), VERSION 0x%08X, STATUS 0x%02X" $id \
          [expr {$id == 0x47534C31 ? "ok" : "WRONG"}] $ver $st]
}

proc gsl_load {file} {
    global GSL_SDRAM GSL_REC_MAGIC
    if {![file exists $file]} { error "$file not found (python3 board_check.py prepare)" }
    set n [expr {[file size $file] / 8}]
    set c [gsl_master]
    master_write_32 $c $GSL_SDRAM 0
    set t0 [clock milliseconds]
    master_write_from_file $c $file [expr {$GSL_SDRAM + 0x100}]
    set ms [expr {[clock milliseconds] - $t0}]
    master_write_32 $c [expr {$GSL_SDRAM + 4}] $n
    master_write_32 $c $GSL_SDRAM $GSL_REC_MAGIC         ;# last: the Nios V starts on this
    close_service master $c
    puts "loaded $n frames from $file in $ms ms"
}

proc gsl_wait_dump {file {timeout_s 600}} {
    global GSL_SDRAM GSL_LOG_MAGIC GSL_LOG_OFF
    set base [expr {$GSL_SDRAM + $GSL_LOG_OFF}]
    set c [gsl_master]
    set t0 [clock seconds]
    while {[gsl_rd $c $base] != $GSL_LOG_MAGIC} {
        if {[clock seconds] - $t0 > $timeout_s} {
            close_service master $c
            error "no result log after $timeout_s s"
        }
        after 250
    }
    set nw [gsl_rd $c [expr {$base + 4}]]
    set nbytes [expr {64 + 4 * $nw}]
    file mkdir [file dirname $file]
    master_read_to_file $c $file $base $nbytes
    master_write_32 $c $base 0
    close_service master $c
    puts "saved $nbytes bytes ($nw words) to $file"
}

proc gsl_run {name} {
    gsl_load rec/mb_$name.bin
    gsl_wait_dump log/$name.bin
    puts "now: python3 board_check.py check $name"
}
