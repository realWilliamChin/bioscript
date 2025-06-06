#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from datetime import datetime
from loguru import logger

def setup_logger(log_name):
    """
    配置日志记录器
    
    参数:
        log_name (str): 日志文件名
    """
    # 创建日志目录
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # 生成日志文件名
    log_file = os.path.join(log_dir, f"{log_name}")
    
    # 移除默认的处理器
    logger.remove()
    
    # 添加控制台输出
    logger.add(sys.stderr, level="ERROR")
    
    # 添加文件输出
    logger.add(log_file, rotation="500 MB", level="DEBUG", encoding="utf-8")
    
    return logger

def get_logger(log_name):
    """
    获取配置好的logger实例
    
    参数:
        log_name (str): 日志文件名
    """
    return setup_logger(log_name) 