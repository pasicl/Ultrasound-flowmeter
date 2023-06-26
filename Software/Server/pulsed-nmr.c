#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#define I2C_SLAVE       0x0703 /* Use this slave address */
#define I2C_SLAVE_FORCE 0x0706 /* Use this slave address, even if it
                                  is already in use by a driver! */

#define ADDR_DAC 0x0F /* AD5622 address 3 */

ssize_t i2c_write_data16(int fd, uint16_t data)
{
  uint8_t buffer[2];
  buffer[0] = data >> 8;
  buffer[1] = data;
  return write(fd, buffer, 2);
}

int main(int argc, char *argv[])
{

  int fd, sock_server, sock_client;
  int i2c_fd, i2c_dac;
  void *cfg, *sts;
  volatile uint32_t *tx_data, *rx_freq, *tx_freq;
  volatile uint16_t *rx_rate, *rx_cntr, *tx_cntr;
  volatile int16_t *tx_level;
  volatile uint8_t *rx_rst, *tx_rst;
  volatile uint64_t *rx_data;
  struct sockaddr_in addr;
  uint64_t command, code, data, phase, level, counter;
  uint32_t *pulses;
  uint64_t *buffer;
  int i, n, position, size, yes = 1;
  int ca=0;
  FILE *f1;
  char name_file[8] = "file1";
 
  size = 0;
  pulses = malloc(16777216);
  buffer = malloc(32768);

  if(pulses == NULL || buffer == NULL)
  {
    perror("malloc");
    return EXIT_FAILURE;
  }

  if((fd = open("/dev/mem", O_RDWR)) < 0)
  {
    perror("open");
    return EXIT_FAILURE;
  }

  i2c_dac = 0;
  if((i2c_fd = open("/dev/i2c-0", O_RDWR)) >= 0)
  {
    if(ioctl(i2c_fd, I2C_SLAVE, ADDR_DAC) >= 0)
    {
      if(i2c_write_data16(i2c_fd, 0x0000) > 0)
      {
        i2c_dac = 1;
      }
    }
  }

  sts = mmap(NULL, sysconf(_SC_PAGESIZE), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0x40000000);
  cfg = mmap(NULL, sysconf(_SC_PAGESIZE), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0x40001000);
  rx_data = mmap(NULL, 16*sysconf(_SC_PAGESIZE), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0x40010000);
  tx_data = mmap(NULL, 16*sysconf(_SC_PAGESIZE), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0x40020000);

  rx_rst = ((uint8_t *)(cfg + 0));
  rx_freq = ((uint32_t *)(cfg + 4));
  rx_rate = ((uint16_t *)(cfg + 8));
  tx_level = ((int16_t *)(cfg + 10));
  rx_cntr = ((uint16_t *)(sts + 12));

  tx_rst = ((uint8_t *)(cfg + 1));
  tx_freq = ((uint32_t *)(cfg + 12));
  tx_cntr = ((uint16_t *)(sts + 14));

  *rx_rst |= 1;
  *rx_rst &= ~2;
  *tx_rst |= 1;

  /* set default RX phase increment */
  *rx_freq = (uint32_t)floor(19000000 / 125.0e6 * (1<<30) + 0.5);
  /* set default RX sample rate */
  *rx_rate = 250;

  /* set default TX level */
  *tx_level = 0;

  /* set default TX phase increment */
  *tx_freq = (uint32_t)floor(19000000 / 125.0e6 * (1<<30) + 0.5);

  if((sock_server = socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    perror("socket");
    return EXIT_FAILURE;
  }

  setsockopt(sock_server, SOL_SOCKET, SO_REUSEADDR, (void *)&yes , sizeof(yes));

  /* setup listening address */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(1001);

  if(bind(sock_server, (struct sockaddr *)&addr, sizeof(addr)) < 0)
  {
    perror("bind");
    return EXIT_FAILURE;
  }

  listen(sock_server, 1024);

  while(1)
  {
    if((sock_client = accept(sock_server, NULL, NULL)) < 0)
    {
      perror("accept");
      return EXIT_FAILURE;
    }
    while(1)
    {
      if(recv(sock_client, (char *)&command, 8, MSG_WAITALL) <= 0) break;
      code = command >> 60;
      data = command & 0xfffffffffffffffULL;
      switch(code)
      {
        case 0:
          /* set RX phase increment */
          if(data > 62500000) continue;
          *rx_freq = (uint32_t)floor(data / 125.0e6 * (1<<30) + 0.5);
          break;
        case 1:
          /* set TX phase increment */
          if(data > 62500000) continue;
          *tx_freq = (uint32_t)floor(data / 125.0e6 * (1<<30) + 0.5);
          break;
        case 2:
          /* set RX sample rate */
          if(data <50  || data > 2500) continue;
          *rx_rate = data;
          break;
        case 3:
          /* set TX level */
          if(data > 32766) continue;
          *tx_level = data;
          break;
        case 4:
          /* set pin */
          if(data < 1 || data > 7) continue;
          *tx_rst |= (1 << data);
          break;
        case 5:
          /* clear pin */
          if(data < 1 || data > 7) continue;
          *tx_rst &= ~(1 << data);
          break;
        case 6:
          /* set DAC */
          if(i2c_dac == 0) continue;
          ioctl(i2c_fd, I2C_SLAVE, ADDR_DAC);
          i2c_write_data16(i2c_fd, data);
          break;
        case 7:
          /* clear pulses */
          size = 0;
          break;
        case 8:
          /* add pulse */
          if(size >= 1048576) continue;
          ++size;
          memset(pulses + (size - 1) * 4, 0, 16);
          /* set pulse width */
          memcpy(pulses + (size - 1) * 4, &data, 8);
	  printf("pulse_length = %llu \n",data); 
          break;
        case 9:
          /* set pulse phase and level */
          phase = data >> 16;
          level = data & 0xffff;
          if(phase < 360) pulses[(size - 1) * 4 + 2] = (uint32_t)floor(phase / 360.0 * (1<<30) + 0.5);
          if(level < 32767) pulses[(size - 1) * 4 + 3] = level;
          break;
        case 10:
          /* start sequence */
          printf("Starting sequence \n");
	  printf("data = %llu \n", data);
	  counter = 0;
          position = 0;
          n = 2048 ;

          /* stop RX and TX */
          *rx_rst &= ~2;

          /* clear RX FIFO */
          *rx_rst |= 1; *rx_rst &= ~1;

          /* clear TX FIFO */
          *tx_rst |= 1; *tx_rst &= ~1;

          while(counter < data)
          {
            /* read I/Q samples from RX FIFO */
            if(n > data - counter) n = data - counter;
            if(*rx_cntr < n * 4) usleep(500);
	    printf("counter = %llu, n = %d, rx_cntr = %d\n", counter, n, *rx_cntr);
            /* read I/Q samples from RX FIFO */
            if(*rx_cntr >= n * 4)
            {
/*	      f1 = fopen(name_file, "w");
	      if(f1 == NULL)
		{printf("File could not be opened");
		exit(-1);
	        }
	      printf("print of transfer");*/
	      for(i = 0; i < n * 2; ++i){
		buffer[i] = *rx_data;
//		fprintf(f1,"%llu ", buffer[i]); 	
		}
//		fprintf(f1, "\n");
//		fclose(f1);
		if(send(sock_client, buffer, n * 16, MSG_NOSIGNAL) < 0){
		printf("send failed!\n");	
		 break;
		}
              counter += n;
            }

            /* write pulses to TX FIFO */
	    printf("size = %d \n", size);
            while(*tx_cntr < 16384 && position < size * 4)
            {
	      printf("Send pulse to TX FIFO, position = %d. tx_cntr = %d \n", position,*tx_cntr);
	      printf("Data sent in position %d is %d \n", position, pulses[position]);
              *tx_data = pulses[position];
              ++position;
            }
		if (*tx_cntr >= 16384)
			printf("warning tx_ctr bigger thna 16384\n");
            /* start RX and TX */
            *rx_rst |= 2;
          }

          *rx_rst |= 1;
          *rx_rst &= ~2;
          *tx_rst |= 1;
          break;
      }
    }

    close(sock_client);
  }

  close(sock_server);

  return EXIT_SUCCESS;
}
