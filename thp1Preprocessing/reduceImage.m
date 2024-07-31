function thisimage=reduceImage(thatImage)
    thisImage=downsample(thatImage',2);
    thisimage=downsample(thisImage',2);