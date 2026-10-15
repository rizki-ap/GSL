/* ============================================================
 * GSL board bring-up (Milestone C3) — Nios V/g, bare metal (HAL)
 *
 *  1. Checks the detector's ID register and tests a scratch area of SDRAM.
 *  2. Waits until the PC loads a recording into SDRAM (System Console,
 *     gsl_load in gsl_sc.tcl), replays it through the detector, and copies
 *     every mailbox word into an SDRAM log.
 *  3. Writes the log header; the PC dumps the log (gsl_wait_dump) and
 *     board_check.py compares it with the simulation.
 *
 *  SDRAM layout (offsets from SDRAM_BASE):
 *    0x0000000  recording header: magic "GSLR", frame count
 *    0x0000100  frames, 8 bytes each (up to ~32 MB = 40 s)
 *    0x2000000  log header (16 words: magic "GSLL", word count, status, ...)
 *    0x2000040  log data: raw mailbox words
 *    0x3F00000  last 1 MB: scratch area for the SDRAM test
 * ============================================================ */
#include <stdio.h>
#include <stdint.h>
#include "system.h"
#include "sys/alt_cache.h"

#define GSL_REG(o)  (*(volatile uint32_t *)(uintptr_t)(GSL_BASE + (o)))
#define GSL_RING(i) (*(volatile uint32_t *)(uintptr_t)(GSL_BASE + 0x8000u + 4u * ((i) & 0x1FFFu)))
#define SD(o)       (*(volatile uint32_t *)(uintptr_t)(SDRAM_BASE + (o)))
#define LEDS        (*(volatile uint32_t *)(uintptr_t)(LED_BASE))

/* gsl_mm registers (c1_interface.md section 3) */
#define R_ID      0x000u
#define R_VERSION 0x004u
#define R_CTRL    0x008u
#define R_STATUS  0x00Cu
#define R_LATE    0x018u
#define R_OVF     0x01Cu
#define R_DROPS   0x020u
#define R_EVCNT   0x024u
#define R_RESCNT  0x028u
#define R_MBRD    0x034u
#define R_MBLEVEL 0x038u
#define R_RPBASE  0x040u
#define R_RPN     0x044u
#define R_RPPOS   0x048u

#define CTRL_RUN    0x01u
#define CTRL_REPLAY 0x02u
#define CTRL_SRST   0x04u
#define ST_DONE     0x40u
#define ST_BUSY     0x20u

#define REC_MAGIC  0x47534C52u          /* "GSLR" */
#define LOG_MAGIC  0x47534C4Cu          /* "GSLL" */
#define FRAMES_OFF 0x00000100u
#define LOG_OFF    0x02000000u
#define LOG_DATA   (LOG_OFF + 64u)
#define TEST_OFF   0x03F00000u
#define TEST_WORDS (0x00100000u / 4u)
#define LOG_MAXW   ((TEST_OFF - LOG_DATA) / 4u)

static unsigned long ul(uint32_t v) { return (unsigned long)v; }

static int sdram_test(void)
{
    uint32_t i, err = 0;
    for (i = 0; i < TEST_WORDS; i++)
        SD(TEST_OFF + 4u * i) = (i * 0x9E3779B1u) ^ 0xA5A5A5A5u;
    alt_dcache_flush_all();
    for (i = 0; i < TEST_WORDS; i++)
        if (SD(TEST_OFF + 4u * i) != ((i * 0x9E3779B1u) ^ 0xA5A5A5A5u))
            err++;
    printf("SDRAM test (1 MB at offset 0x%08lx): %s (%lu errors)\n",
           ul(TEST_OFF), err ? "FAILED" : "ok", ul(err));
    return err == 0;
}

static void run_replay(uint32_t nframes)
{
    uint32_t rd = 0, nw = 0, idle = 0, lvl, k, st;

    /* stop, soft reset (clears detector state, mailbox, counters) */
    GSL_REG(R_CTRL) = CTRL_REPLAY | CTRL_SRST;
    while (GSL_REG(R_CTRL) & CTRL_SRST)
        ;
    SD(LOG_OFF) = 0;
    GSL_REG(R_RPBASE) = (uint32_t)SDRAM_BASE + FRAMES_OFF;
    GSL_REG(R_RPN)    = nframes;
    GSL_REG(R_MBRD)   = 0;
    GSL_REG(R_CTRL)   = CTRL_RUN | CTRL_REPLAY;
    LEDS = 0x02;

    for (;;) {
        lvl = GSL_REG(R_MBLEVEL);
        if (lvl) {
            for (k = 0; k < lvl; k++) {
                uint32_t w = GSL_RING(rd + k);
                if (nw < LOG_MAXW)
                    SD(LOG_DATA + 4u * nw) = w;
                nw++;
            }
            rd = (rd + lvl) & 0x1FFFu;
            GSL_REG(R_MBRD) = rd;
            idle = 0;
        } else {
            st = GSL_REG(R_STATUS);
            if ((st & ST_DONE) && !(st & ST_BUSY)) {
                if (++idle > 2000000u)        /* quiet for ~0.5 s: finished */
                    break;
            } else
                idle = 0;
        }
    }
    GSL_REG(R_CTRL) = CTRL_REPLAY;            /* stop */

    st = GSL_REG(R_STATUS);
    SD(LOG_OFF + 4)  = (nw < LOG_MAXW) ? nw : LOG_MAXW;
    SD(LOG_OFF + 8)  = st;
    SD(LOG_OFF + 12) = GSL_REG(R_LATE);
    SD(LOG_OFF + 16) = GSL_REG(R_OVF);
    SD(LOG_OFF + 20) = GSL_REG(R_DROPS);
    SD(LOG_OFF + 24) = GSL_REG(R_EVCNT);
    SD(LOG_OFF + 28) = GSL_REG(R_RESCNT);
    SD(LOG_OFF + 32) = GSL_REG(R_RPPOS);
    SD(LOG_OFF + 36) = nframes;
    for (k = 40; k < 64; k += 4)
        SD(LOG_OFF + k) = 0;
    alt_dcache_flush_all();
    SD(LOG_OFF) = LOG_MAGIC;                  /* written last: the PC waits for it */
    alt_dcache_flush_all();
    LEDS = 0x04;

    printf("Done: %lu words logged%s. SW %lu, MB %lu, results %lu; late %lu, overflow %lu, "
           "mailbox drops %lu; STATUS 0x%02lx (underflow %lu, result error %lu)\n",
           ul(nw), (nw > LOG_MAXW) ? " (LOG FULL)" : "",
           ul(GSL_REG(R_EVCNT) & 0xFFFFu), ul(GSL_REG(R_EVCNT) >> 16), ul(GSL_REG(R_RESCNT)),
           ul(GSL_REG(R_LATE)), ul(GSL_REG(R_OVF)), ul(GSL_REG(R_DROPS)),
           ul(st), ul(st & 1u), ul((st >> 4) & 1u));
}

int main(void)
{
    uint32_t id, n, blink = 0;

    printf("\nGSL board bring-up (Milestone C3)\n");
    id = GSL_REG(R_ID);
    printf("gsl_mm ID 0x%08lx (%s), VERSION 0x%08lx\n", ul(id),
           id == 0x47534C31u ? "ok" : "WRONG", ul(GSL_REG(R_VERSION)));
    sdram_test();

    for (;;) {
        printf("Waiting for a recording (System Console: gsl_load <file.bin>) ...\n");
        do {
            alt_dcache_flush_all();             /* see the PC's writes to SDRAM */
            LEDS = ((++blink >> 18) & 1u);
        } while (SD(0) != REC_MAGIC);
        n = SD(4);
        printf("Recording: %lu frames (%lu ms). Replaying in real time ...\n",
               ul(n), ul(n / 100u));
        run_replay(n);
        SD(0) = 0;                              /* consumed */
        alt_dcache_flush_all();
    }
    return 0;
}
