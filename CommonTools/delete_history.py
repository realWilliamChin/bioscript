#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2025/12/05 11:48
# Author        : William GoGo
"""
递归删除指定目录下所有 .history 和 Rplots.pdf 文件
"""
import os
import sys
import argparse
from pathlib import Path
from loguru import logger

DEFAULT_TARGET_FILES = ['.history', 'Rplots.pdf']


def should_delete_file(file_path, target_files):
    """判断文件是否应该被删除"""
    # 检查文件名是否在目标文件列表中
    # 或者扩展名为 .history（保持向后兼容）
    return file_path.name in target_files or file_path.suffix == '.history'


def scan_files_to_delete(target_dir, target_files):
    """扫描目录，返回所有需要删除的文件列表"""
    target_path = Path(target_dir)
    
    if not target_path.exists() or not target_path.is_dir():
        logger.error(f"目标目录不存在或不是目录: {target_dir}")
        return []
    
    files_to_delete = []
    try:
        for root, dirs, files in os.walk(target_path):
            for file_name in files:
                file_path = Path(root) / file_name
                if should_delete_file(file_path, target_files):
                    files_to_delete.append(file_path)
    except Exception as e:
        logger.error(f"遍历目录时出错: {e}")
    
    return files_to_delete


def delete_files(files_to_delete, dry_run=False, verbose=False):
    """删除指定的文件列表"""
    deleted_count = error_count = 0
    
    if dry_run:
        logger.warning("*** 模拟运行模式，不会实际删除文件 ***")
    
    for file_path in files_to_delete:
        try:
            if verbose or dry_run:
                logger.info(f"{'[模拟] ' if dry_run else ''}删除: {file_path}")
            if not dry_run:
                file_path.unlink()
            deleted_count += 1
        except (PermissionError, OSError) as e:
            logger.error(f"删除失败 {file_path}: {e}")
            error_count += 1
        except KeyboardInterrupt:
            logger.warning("用户中断操作")
            break
    
    return deleted_count, error_count


def main():
    parser = argparse.ArgumentParser(description='递归删除指定目录下指定的文件（默认: .history 和 Rplots.pdf）')
    parser.add_argument('directory', help='要清理的目标目录路径')
    parser.add_argument('--files', '-f', nargs='+', default=DEFAULT_TARGET_FILES,
                        help=f'要删除的文件名列表（默认: {" ".join(DEFAULT_TARGET_FILES)}）')
    parser.add_argument('--dry-run', action='store_true', help='模拟运行模式，不实际删除')
    parser.add_argument('--verbose', '-v', action='store_true', help='显示详细信息')
    parser.add_argument('--force', action='store_true', help='跳过确认提示')
    
    args = parser.parse_args()
    
    # 将目标文件列表转换为集合以便快速查找
    target_files = set(args.files)
    
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if args.verbose else "INFO")
    
    logger.info(f"开始扫描目录: {Path(args.directory).resolve()}")
    logger.info(f"目标文件: {', '.join(sorted(target_files))}")
    files_to_delete = scan_files_to_delete(args.directory, target_files)
    
    if not files_to_delete:
        logger.info("未找到需要删除的文件")
        return 0
    
    # 打印所有要删除的文件
    logger.info(f"\n找到 {len(files_to_delete)} 个文件需要删除:")
    for file_path in files_to_delete:
        logger.info(f"  {file_path}")
    
    # 如果没有使用 --force，询问用户确认
    if not args.force:
        logger.info(f"\n模式: {'模拟运行' if args.dry_run else '实际删除'}")
        try:
            response = input("\n是否继续删除这些文件? [y/N]: ").strip().lower()
            if response not in ['y', 'yes']:
                logger.info("操作已取消")
                return 0
        except (KeyboardInterrupt, EOFError):
            logger.warning("\n操作已取消")
            return 0
    else:
        logger.info(f"模式: {'模拟运行' if args.dry_run else '实际删除'} (已使用 --force，跳过确认)")
    
    # 执行删除
    deleted_count, error_count = delete_files(files_to_delete, args.dry_run, args.verbose)
    logger.info(f"清理完成 | 删除: {deleted_count} | 错误: {error_count}")
    
    return 0 if error_count == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
