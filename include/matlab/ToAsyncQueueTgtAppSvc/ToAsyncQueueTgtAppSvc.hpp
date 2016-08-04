/* Copyright 2013-2015 MathWorks, Inc. */

#ifndef ToAsyncQueueTgtAppSvc_hpp
#define ToAsyncQueueTgtAppSvc_hpp

#include "ToAsyncQueueTgtAppSvc_dll.hpp"
#include "coder/target_services/Application.hpp"

#ifdef BUILDING_LIBMWCODER_TOASYNCQUEUETGTAPPSVC
#  include "coder/target_services/CommService.hpp"
#else
#  include "CommService.hpp"
#endif

class TOASYNCQUEUETGTAPPSVC_API ToAsyncQueueTgtAppSvc : public coder::tgtsvc::Application
{
  public:
    ToAsyncQueueTgtAppSvc();
    ~ToAsyncQueueTgtAppSvc();

    void sendData(uint32_t id, double time, void *data, uint32_t sizeOfData);
    void handleMessage(coder::tgtsvc::Message *message);

    uint8_t id() { return(coder::tgtsvc::Application::TO_ASYNC_QUEUE_ID); }

    virtual void handleConnect(bool connected) {};

  private:
    ToAsyncQueueTgtAppSvc(const ToAsyncQueueTgtAppSvc &);                 
    const ToAsyncQueueTgtAppSvc& operator=(const ToAsyncQueueTgtAppSvc &);
};

#endif

