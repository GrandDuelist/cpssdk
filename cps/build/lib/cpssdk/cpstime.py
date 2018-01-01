from datetime import datetime
from datetime import timedelta

class CPSTime():
    def parseTime(self,timeStr):
        return datetime.strptime(timeStr,'%Y-%m-%d %H:%M:%S')
    
    def parseTimeToStr(self,time_value):
        return(time_value.strftime("%Y-%m-%d %H:%M:%S"))
    
    def oneDayToSeconds(self,time_str):
        time = datetime.strptime(time_str,'%H:%M:%S')
        time_zero = datetime.strptime('00:00:00','%H:%M:%S')
        return((time-time_zero).total_seconds())
    
    def timeSlot(self,time=None,t_hour = None, t_min=None,t_sec = None):
        divide_min = self.timeToMin(t_hour,t_min,t_sec)
        start_min = self.timeToMin(time.hour,time.minute,time.second)
        slot = int(start_min/divide_min)
        return(slot)
    
    def timeToMin(self,t_hour=None,t_min=None,t_sec=None):
        total_min = 0
        if t_hour is not None:
            total_min = total_min + t_hour * 60
        if t_min is not None:
            total_min = total_min + t_min
        if t_sec is not None:
            total_min = float(total_min) + float(t_sec)/float(60)
        return(total_min)
    
    def timeDifferenceToSeconds(self,start_time,end_time):
        return((end_time - start_time).total_seconds())
    
    def changeTimeSlot(self,time_slot,orig_time_granularity,new_time_granularity):
        new_time_slot = (time_slot * orig_time_granularity)/new_time_granularity
        return(int(new_time_slot))
    
    def timeSlotToTime(self,time_slot,time_slot_granularity=1,start_time = None):
        if start_time is None:
            time_zero = datetime.strptime('00:00:00','%H:%M:%S')
        else:
            time_zero = start_time
        time_dif = timedelta(seconds = time_slot * time_slot_granularity)
        time = time_zero+time_dif
        # return(datetime.time(time))
        return(time)
    
    def equalTimeSeparator(self,start_time,end_time,interval_seconds = 600):
        result = []
        time_dif = timedelta(seconds = interval_seconds)
        one_time = start_time
        while(one_time < end_time):
            result.append(one_time)
            one_time = one_time + time_dif
        result.append(one_time)
        return(result)
