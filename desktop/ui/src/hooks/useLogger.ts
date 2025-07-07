import { useCallback } from 'react'

interface LogLevel {
  DEBUG: 'debug'
  INFO: 'info'
  WARN: 'warn'
  ERROR: 'error'
}

const LOG_LEVELS: LogLevel = {
  DEBUG: 'debug',
  INFO: 'info',
  WARN: 'warn',
  ERROR: 'error'
}

export const useLogger = () => {
  const log = useCallback((level: keyof LogLevel, message: string, data?: any) => {
    const timestamp = new Date().toISOString()
    const logMessage = `[${timestamp}] [${level}] ${message}`
    
    switch (level) {
      case 'DEBUG':
        console.debug(logMessage, data)
        break
      case 'INFO':
        console.info(logMessage, data)
        break
      case 'WARN':
        console.warn(logMessage, data)
        break
      case 'ERROR':
        console.error(logMessage, data)
        break
      default:
        console.log(logMessage, data)
    }
  }, [])

  const debug = useCallback((message: string, data?: any) => {
    log('DEBUG', message, data)
  }, [log])

  const info = useCallback((message: string, data?: any) => {
    log('INFO', message, data)
  }, [log])

  const warn = useCallback((message: string, data?: any) => {
    log('WARN', message, data)
  }, [log])

  const error = useCallback((message: string, data?: any) => {
    log('ERROR', message, data)
  }, [log])

  return {
    log,
    debug,
    info,
    warn,
    error
  }
} 