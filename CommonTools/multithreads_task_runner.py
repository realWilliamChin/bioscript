#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2026/03/16 13:30
# Author        : William GoGo
import concurrent.futures
from loguru import logger
from typing import List, Callable, Any, Tuple


def run_multithreads_tasks(
    task_func: Callable,
    tasks_args: List[Tuple],
    max_workers: int = 3,
    show_log: bool = True
) -> List[Any]:
    """
    通用多进程任务运行器，保持固定数量的并行任务，动态提交新任务
    
    参数:
        task_func: 需要执行的任务函数
        tasks_args: 所有任务的参数列表，每个元素是一个 tuple，包含对应任务函数的所有参数
        max_workers: 最大并行任务数，默认 3
        show_log: 是否显示运行日志，默认 True
    
    返回:
        results: 所有任务的返回结果列表，顺序和任务参数列表顺序一致
    """
    results = [None] * len(tasks_args)
    
    # 控制同时运行的任务数量
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 活跃的任务列表，存储 (future, task_index)
        active_futures = []
        # 任务索引
        task_index = 0
        total_tasks = len(tasks_args)
        
        # 首先提交 max_workers 个任务
        while len(active_futures) < max_workers and task_index < total_tasks:
            args = tasks_args[task_index]
            if show_log:
                logger.info(f'开始执行第 {task_index + 1}/{total_tasks} 个任务，参数: {args}')
            
            future = executor.submit(task_func, *args)
            active_futures.append((future, task_index))
            task_index += 1
        
        # 处理完成的任务并提交新任务
        while active_futures:
            # 等待任意一个任务完成
            done, not_done = concurrent.futures.wait(
                [f for f, _ in active_futures],
                return_when=concurrent.futures.FIRST_COMPLETED
            )
            
            # 处理完成的任务
            for future in done:
                # 找到对应的任务索引
                for i, (f, idx) in enumerate(active_futures):
                    if f == future:
                        task_index_original = idx
                        active_futures.pop(i)
                        break
                
                try:
                    result = future.result()
                    results[task_index_original] = result
                    if show_log:
                        logger.info(f'第 {task_index_original + 1} 个任务执行完成，结果: {result}')
                except Exception as exc:
                    if show_log:
                        logger.error(f'第 {task_index_original + 1} 个任务运行时出错: {exc}')
                    results[task_index_original] = None
                
                # 如果还有未处理的任务，提交新任务
                if task_index < total_tasks:
                    new_args = tasks_args[task_index]
                    if show_log:
                        logger.info(f'开始执行第 {task_index + 1}/{total_tasks} 个任务，参数: {new_args}')
                    
                    new_future = executor.submit(task_func, *new_args)
                    active_futures.append((new_future, task_index))
                    task_index += 1
    
    return results


def run_multithreads_threads(
    task_func: Callable,
    tasks_args: List[Tuple],
    max_workers: int = 3,
    show_log: bool = True
) -> List[Any]:
    """
    通用多线程任务运行器，保持固定数量的并行任务，动态提交新任务
    使用 ThreadPoolExecutor 适用于 I/O 密集型任务
    
    参数:
        task_func: 需要执行的任务函数
        tasks_args: 所有任务的参数列表，每个元素是一个 tuple，包含对应任务函数的所有参数
        max_workers: 最大并行任务数，默认 3
        show_log: 是否显示运行日志，默认 True
    
    返回:
        results: 所有任务的返回结果列表，顺序和任务参数列表顺序一致
    """
    results = [None] * len(tasks_args)
    
    # 控制同时运行的任务数量
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 活跃的任务列表，存储 (future, task_index)
        active_futures = []
        # 任务索引
        task_index = 0
        total_tasks = len(tasks_args)
        
        # 首先提交 max_workers 个任务
        while len(active_futures) < max_workers and task_index < total_tasks:
            args = tasks_args[task_index]
            if show_log:
                logger.info(f'开始执行第 {task_index + 1}/{total_tasks} 个任务，参数: {args}')
            
            future = executor.submit(task_func, *args)
            active_futures.append((future, task_index))
            task_index += 1
        
        # 处理完成的任务并提交新任务
        while active_futures:
            # 等待任意一个任务完成
            done, not_done = concurrent.futures.wait(
                [f for f, _ in active_futures],
                return_when=concurrent.futures.FIRST_COMPLETED
            )
            
            # 处理完成的任务
            for future in done:
                # 找到对应的任务索引
                for i, (f, idx) in enumerate(active_futures):
                    if f == future:
                        task_index_original = idx
                        active_futures.pop(i)
                        break
                
                try:
                    result = future.result()
                    results[task_index_original] = result
                    if show_log:
                        logger.info(f'第 {task_index_original + 1} 个任务执行完成，结果: {result}')
                except Exception as exc:
                    if show_log:
                        logger.error(f'第 {task_index_original + 1} 个任务运行时出错: {exc}')
                    results[task_index_original] = None
                
                # 如果还有未处理的任务，提交新任务
                if task_index < total_tasks:
                    new_args = tasks_args[task_index]
                    if show_log:
                        logger.info(f'开始执行第 {task_index + 1}/{total_tasks} 个任务，参数: {new_args}')
                    
                    new_future = executor.submit(task_func, *new_args)
                    active_futures.append((new_future, task_index))
                    task_index += 1
    
    return results