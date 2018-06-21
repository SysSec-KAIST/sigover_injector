typedef struct {
  char *devname;
  uhd_usrp_handle usrp;
  uhd_rx_streamer_handle rx_stream;
  uhd_tx_streamer_handle tx_stream;

  uhd_rx_metadata_handle rx_md, rx_md_first;
  uhd_tx_metadata_handle tx_md;

  uhd_meta_range_handle rx_gain_range;
  size_t rx_nof_samples;
  size_t tx_nof_samples;
  double tx_rate;
  bool dynamic_rate;
  bool has_rssi;
  uint32_t nof_rx_channels;
  int nof_tx_channels;

  srslte_rf_error_handler_t uhd_error_handler;

  float current_master_clock;

  bool async_thread_running;
  pthread_t async_thread;
} rf_uhd_handler_t;
