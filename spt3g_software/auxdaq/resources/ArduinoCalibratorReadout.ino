//#define DEBUG_SERIAL
#define INCLUDE_ERRORS

#include <SPI.h>
#include <Ethernet2.h>
#include <EthernetUdp2.h>   

#define MAX_LOOP_TIME 0xC000

IPAddress ip(192, 168, 1, 177); //ip address of the Arduino
IPAddress remote_ip(192, 168, 1, 42); //ip address we are sending it to

byte mac[] = { 0x90, 0xA2, 0xDA, 0x10, 0x3B, 0x16 };
unsigned int localPort = 42668;
EthernetUDP Udp;

/**
Generic Notes:


No double buffering

Timer 3 is the clock
Timer 1 is used to trigger checking calibrator by just haveing it count forever and trigger an interrupt on overflow
**/


/**
Constants selecting hardware configuration
**/


//pin 13 is the input capture for counter 3 which we use for the cpu clock, (feed in the irig information)
//pin 0 is used as the input for the calibrator

// So, because of insanity, Pin 0 as labelled on the arduino is PortD bit 2 (0 indexed) on the atmega processor.
//  This means all of the interrupts and software readout uses bit 2, while it's plugged in to pin 0.


/**
Infor we will be sending back
**/

struct
CalibratorInfo{
  uint64_t clock_cnt; //clock count
  unsigned int is_high;
};

struct IrigInfo{
  unsigned int random_header = 0xCAFE; //I need my coffee
  uint64_t rising_edge_time; //
  unsigned int info[10];
  uint64_t re_count[10];  
};

/**
Interrupts for controlling the over flow of timer 3 whcih is handling the clock
**/

volatile uint64_t counter_3_overflow = 0;

//sets the overflow procedure for timer 3
ISR(TIMER3_OVF_vect){
  counter_3_overflow++;
}

/**
Error reporting
**/
#define ERR_NONE 0
#define ERR_DESYNC 1
#define ERR_IRIG_NO_SEND 2
#define ERR_IRIG_BIT_WEIRD 4
#define ERR_BAD_INTERRUPT 8           
#define ERR_GLITCH 16

struct ErrorInfo {
  unsigned int header = 0xE12A;
  unsigned int err_code = ERR_NONE;  
};

volatile struct ErrorInfo error_state;

/**
Generic IRIG Routines
**/

//Define the possible PWM irig bits, 0,1, PI is the synchronization pulse, ERR is an error...
#define IRIG_0 0
#define IRIG_1 1
#define IRIG_PI 2
#define IRIG_ERR 3
#define IRIG_GLITCH 4

struct IrigInfo irig_packet_0;
struct IrigInfo irig_packet_1;

volatile struct IrigInfo * the_irig_packet = &irig_packet_0;
volatile struct IrigInfo * to_send_irig_packet = &irig_packet_1;
volatile bool irig_packet_pointer_coord = false;

/**
Variables for Irig Interrupt
**/

uint64_t rising_edge_t;
uint64_t falling_edge_t;

volatile unsigned char prev_bit_type = IRIG_ERR;
volatile unsigned char irig_parser_is_synched = 0;
volatile unsigned char bit_position = 0;

volatile unsigned char book_the_irig_packet = 0;

void irig_interrupt(){
  //Grab the clock
  unsigned char overflow_bits = TIFR3;
  unsigned int clock_cnt = ICR3;
  uint64_t counter_3_overflow_cached = counter_3_overflow;
  
  sei(); //enable other itnerrupts so we sample the encoder regularly
  
  overflow_bits = (overflow_bits & 1) && (clock_cnt < MAX_LOOP_TIME);
  unsigned char irig_bit_type = IRIG_ERR;
  
  TCCR3B ^= 1 << 6; //toggle the interrupt on the capture flag for rising or falling
  
  
  //if falling edge and rising edge valid, interpret bit type
  //if (is_falling_edge){
  if (TCCR3B & 64) {
     #ifdef DEBUG_SERIAL
     Serial.println("At falling edge:");
     #endif
     falling_edge_t = clock_cnt + (((uint64_t)(overflow_bits + counter_3_overflow_cached)) << 16);

     uint64_t delta = falling_edge_t - rising_edge_t; 
     if (delta < 1000) irig_bit_type = IRIG_GLITCH;
     else if (delta < 56000) irig_bit_type = IRIG_0;
     else if ( delta > 196608) irig_bit_type = IRIG_ERR;
     else if ( delta > 104000 ) irig_bit_type = IRIG_PI;
     else irig_bit_type = IRIG_1;

     if (irig_bit_type == IRIG_GLITCH) {
        error_state.err_code = ERR_GLITCH;
        return; 
     }

     if (irig_bit_type == IRIG_ERR) {
           irig_parser_is_synched = 0;
           error_state.err_code = ERR_IRIG_BIT_WEIRD;
           bit_position = 0;
     }
     // if we have synched actually parse the information from this bit
     if (irig_parser_is_synched) {
       #ifdef DEBUG_SERIAL
       Serial.println("synched falling edge");
       #endif   
       if ( bit_position == 100 ) {
         irig_packet_pointer_coord = !irig_packet_pointer_coord;
         if (irig_packet_pointer_coord) {
           the_irig_packet = &irig_packet_1;
           to_send_irig_packet = &irig_packet_0;
         } else {
           the_irig_packet = &irig_packet_0;
           to_send_irig_packet = &irig_packet_1;
         }      
         book_the_irig_packet = 1;
         bit_position = 1;
         the_irig_packet->rising_edge_time = rising_edge_t;
       } else if ( bit_position % 10 == 9) {
         if (irig_bit_type != IRIG_PI) {
           bit_position = 0;
           irig_parser_is_synched = 0;
           error_state.err_code = ERR_DESYNC;
         }
         unsigned char ind = bit_position/10;
         the_irig_packet->re_count[ind] = rising_edge_t;
         bit_position++;
       } else if (irig_bit_type != IRIG_0 && irig_bit_type != IRIG_1) {
             bit_position = 0;
             irig_parser_is_synched = 0;
             error_state.err_code = ERR_DESYNC | ERR_IRIG_BIT_WEIRD;
        } else {
         unsigned char offset = bit_position % 10;
         the_irig_packet->info[bit_position/10] &= ~(1 << offset);
         the_irig_packet->info[bit_position/10] |= irig_bit_type << (offset);   
         bit_position++;
       }   
     } else {
       #ifdef DEBUG_SERIAL
       Serial.println("not synched falling edge");
       #endif
       //if we aren't synchronized we look for the start of an irig frame.
       if ( irig_bit_type == IRIG_PI && prev_bit_type == IRIG_PI) {
         bit_position = 1;
         irig_parser_is_synched = 1;
         the_irig_packet->rising_edge_time = rising_edge_t;
       }
     }
     prev_bit_type = irig_bit_type;
   } else {
     rising_edge_t = clock_cnt + (((uint64_t)(overflow_bits + counter_3_overflow_cached)) << 16);
   } 
}

ISR(TIMER3_CAPT_vect){
  irig_interrupt();
}


#define NCAL 4
#define HALF_NCAL 2

unsigned int calibrator_info_header = 0xBBBB;

volatile unsigned char send_calibrator_data = 0;
volatile unsigned char send_cal_start_ind = 0;

struct CalibratorInfo cinfo[NCAL];
unsigned char curr_cal = 0;

void calibrator_interrupt(){
  cinfo[curr_cal].clock_cnt = TCNT3;
  unsigned char overflow_bits = TIFR3;
  cinfo[curr_cal].clock_cnt +=  (counter_3_overflow + 
       ((overflow_bits & 1) && (cinfo[curr_cal].clock_cnt < MAX_LOOP_TIME))) << 16;
  cinfo[curr_cal].is_high = (PIND & 4) >> 2;
  
  
  curr_cal = (curr_cal+1) % NCAL;
  
  if (curr_cal % HALF_NCAL == 0) {
    send_calibrator_data = 1;
    send_cal_start_ind = curr_cal == 0 ? HALF_NCAL : 0;
  }  
}  

ISR(INT2_vect){  
  calibrator_interrupt(); 
}


ISR(TIMER1_OVF_vect){ 
  calibrator_interrupt(); 
}

void setup() {    
  Ethernet.begin(mac, ip);
  Udp.begin(localPort);
  // register the irig terminal as an input
  TCCR1A = 0; //
  TCCR1B = 2; //tells timer 1 to use clk 1 as a decimated input source
  TIMSK1 = 1; //Sets counter 1 to set interrupt on overflow
 
  TCCR3A = 0; //
  TCCR3B = 65; //set counter 3 to use the clock as it's source and use rising edge to trigger input capture
  TIMSK3 = 33; //Sets counter 1 to send an interrupt on overflow and interrupt on input capture
  
  //set up calibrator input
  DDRD = 0;
  PORTD = 0;

  EICRA  &= ~(1 << 5); //sets it so the logic change causes an interrupt
  EICRA  |= 1 << 4;   
  EIMSK |= 1 << 2;//sets interrupt enabled for the calibrator pin
  
  for (int i=0; i < 10; i++) the_irig_packet->info[i] = 0;
}





void loop() {

  if (book_the_irig_packet) {
    Udp.beginPacket(remote_ip, localPort);
    Udp.write(((char*)to_send_irig_packet), sizeof(struct IrigInfo));
    Udp.endPacket();
    book_the_irig_packet = 0;
  }
  if (send_calibrator_data){
    Udp.beginPacket(remote_ip, localPort);
    Udp.write((char*)&calibrator_info_header, sizeof(unsigned int));
    Udp.write((char*)(cinfo + send_cal_start_ind), HALF_NCAL * sizeof(CalibratorInfo));
    Udp.endPacket();
    send_calibrator_data = 0; 
  }
  
  if (error_state.err_code) {
    Udp.beginPacket(remote_ip, localPort);
    Udp.write(((char*)&error_state), sizeof(struct ErrorInfo));
    Udp.endPacket();  
    error_state.err_code = ERR_NONE;
  }
}



